# author: "Angie Madrid"
# date: "December 17, 2015"

#-------------------------
# libraries
#-------------------------
library(caret)
library(e1071)
library(pROC)

#-------------------------
# LOAD DATA AND PREPROCESS
#-------------------------
data.raw<-read.csv("alldata_genes.csv")
summary(data.raw)
data.raw$seg_mean<-round(data.raw$seg_mean)
seg_mean.max<-max(data.raw$seg_mean)
seg_mean.min<-min(data.raw$seg_mean)
data.raw$cancer<-factor(data.raw$cancer,levels=c("1","-1"))
data.raw$seg_size<-abs(data.raw$start-data.raw$end)

#create binary matrix of chromosome data and concat with other variables 
data<-with(data.raw,data.frame(model.matrix(~chromosome-1,data.raw),
                               num_probes,seg_mean,num_genes,
                               most_freq,seg_size,cancer))

# FUNCTION
# normalize the data
normalize<- function(x) {
  return((x-min(x))/(max(x)-min(x)))
}

data$num_probes<-normalize(data$num_probes)
data$num_genes<-normalize(data$num_genes)
data$most_freq<-normalize(data$most_freq)
data$seg_size<-normalize(data$seg_size)
levels(data$cancer)<-c("cancer","no_cancer")

#-------------------------
# SETUP TRAINING AND TESTING DATA SETS
#-------------------------
data.len<-sample(1:nrow(data),3*nrow(data)/4)
data.training<-data[data.len,]
data.testing<-data[-data.len,]
length(which(data.training$cancer=="cancer"))
length(which(data.training$cancer=="no_cancer"))
length(which(data.testing$cancer=="cancer"))
length(which(data.testing$cancer=="no_cancer"))
write.csv(data.training,"training.csv",row.names=FALSE,col.names = FALSE)
write.csv(data.testing,"testing.csv",row.names=FALSE,col.names = FALSE)

fitControl <- trainControl(method = "repeatedcv",
                           number = 5,
                           repeats = 5,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

#-------------------------
# BUILD THE RANDOM FOREST MODEL
#-------------------------
# 24 hours to run
rfFit <- train(cancer ~ ., data = data.training,
               method = "rf",
               trControl = fitControl,
               ## Specify which metric to optimize
               metric = "ROC")
#importance
rfImp<-varImp(rfFit,scale=FALSE)
rfImp
#  only 20 most important variables shown (out of 29)
#
#             Overall
#seg_size     4421.59
#num_probes   3071.26
#num_genes    2014.28
#most_freq    1069.11
#seg_mean      478.61
#chromosome1   146.92
#chromosome6   146.69
#chromosome8   143.94
#chromosome2   140.85
#chromosome11  138.18
#chromosome5   136.22
#chromosome3   128.06
#chromosome10  124.92
#chromosome12  114.94
#chromosomeX   111.51
#chromosome9   110.95
#chromosome4   110.47
#chromosome16  108.60
#chromosome14  102.76
#chromosome13   99.43

plot(rfFit,metric="ROC",main="Random Forest:5-fold Cross Validation")

#-------------------------
# PREDICT USING THE RF MODEL
#-------------------------
rf.pred<-predict(rfFit,data.testing)
rf.ctable<-confusionMatrix(rf.pred,data.testing$cancer)
rf.ctable
#Confusion Matrix and Statistics
#
#           Reference
#Prediction  cancer no_cancer
#  cancer     11976      2363
#  no_cancer    848      1076
#                                          
#               Accuracy : 0.8026          
#                 95% CI : (0.7964, 0.8087)
#    No Information Rate : 0.7885          
#    P-Value [Acc > NIR] : 5.332e-06       
#                                          
#                  Kappa : 0.2942          
# Mcnemar's Test P-Value : < 2.2e-16       
#                                          
#            Sensitivity : 0.9339          
#            Specificity : 0.3129          
#         Pos Pred Value : 0.8352          
#         Neg Pred Value : 0.5593          
#             Prevalence : 0.7885          
#         Detection Rate : 0.7364          
#   Detection Prevalence : 0.8817          
#      Balanced Accuracy : 0.6234          
#                                          
#       'Positive' Class : cancer          

rf.acc<-rf.ctable$overall
rf.acc
#      Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#  8.025580e-01   2.941765e-01   7.963554e-01   8.086521e-01   7.885384e-01   5.331754e-06  2.899818e-157 

rfROC<-roc(data.testing$cancer, as.numeric(rf.pred),levels=c("cancer","no_cancer"))
plot(rfROC,main="Random Forest: 5-fold cross validation")
auc(rfROC)
#Area under the curve: 0.6234


#-------------------------
# BUILD THE RANDOM FOREST MODEL
# NOT USING CARET
# https://chandramanitiwary.wordpress.com/2014/03/17/r-tips-part2-rocr-example-with-randomforest/
#-------------------------
library(ROCR)
rf.bagging<-randomForest(cancer~.,data=data.training,importance=TRUE)
rf.bag.acc<-predict(rf.bagging,data.testing)
rf.bag.ctable<-confusionMatrix(rf.bag.acc,data.testing$cancer)
rf.bag.ctable
#Confusion Matrix and Statistics
#
#           Reference
#Prediction  cancer no_cancer
#  cancer     12684      2763
#  no_cancer    140       676
#                                          
#               Accuracy : 0.8215          
#                 95% CI : (0.8155, 0.8274)
#    No Information Rate : 0.7885          
#    P-Value [Acc > NIR] : < 2.2e-16       
#                                          
#                  Kappa : 0.2575          
# Mcnemar's Test P-Value : < 2.2e-16       
#                                          
#            Sensitivity : 0.9891          
#            Specificity : 0.1966          
#         Pos Pred Value : 0.8211          
#         Neg Pred Value : 0.8284          
#             Prevalence : 0.7885          
#         Detection Rate : 0.7799          
#   Detection Prevalence : 0.9498          
#      Balanced Accuracy : 0.5928          
#                                          
#       'Positive' Class : cancer          

#-------------------------
# PREDICT THE RANDOM FOREST MODEL
# NOT USING CARET
#-------------------------
rf.bag.pred<-predict(rf.bagging,data.testing,type='prob')
rf.bag.import<-importance(rf.bagging)
rf.bag.import_reordered <- rf.bag.import[order(rf.bag.import[,4], decreasing=TRUE),]
rf.bag.import_reordered
#                cancer   no_cancer MeanDecreaseAccuracy MeanDecreaseGini
#seg_size     29.961188  32.7086146            49.858443       1705.44001
#num_probes   37.343564  -2.0173547            46.654014       1489.37952
#num_genes    39.502861   0.6679628            48.291341       1094.95600
#most_freq    43.537433 -18.9406917            46.200470        576.30958
#seg_mean     52.710636  -9.7997878            63.619187        333.49818
#chromosome6  26.238913   6.6659771            30.901876         66.84105
#chromosome8  15.335383   5.9579444            19.268143         64.44571
#chromosome11  9.775470   8.2298696            14.256051         61.75437
#chromosome2  21.865717   7.7852181            24.903234         61.28594
#chromosome5  17.446276   7.9920999            20.146874         60.76391
#chromosome12 17.365170  14.0047701            25.531003         57.38079
#chromosome10  6.320144  16.0444302            18.358655         56.75473
#chromosome3  21.397899  14.7040631            27.108491         54.57633
#chromosomeX  23.993604  -7.5897551            19.858749         54.53954
#chromosome1  34.469924  27.6489406            43.247384         53.96470
#chromosome17 16.896747   3.5012195            20.055055         53.14371
#chromosome22 11.687916  27.6951613            30.329016         53.02085
#chromosome14  9.650905  14.0415123            17.083639         51.25456
#chromosome9  15.120663   1.6971299            16.358067         51.02333
#chromosomeY  24.047586  24.4850260            27.953692         49.85505
#chromosome13  9.785086  19.5657320            23.090953         49.12060
#chromosome16 22.858393  11.2921613            28.101172         48.74237
#chromosome4  31.387896  20.3909525            39.098789         46.56755
#chromosome7  30.256888  20.6370656            36.208008         45.86445
#chromosome19 17.875946   1.2834978            19.734361         43.69314
#chromosome20  2.421472   8.4716120             8.115813         42.51536
#chromosome21  7.109922  17.6107216            19.995132         41.96372
#chromosome18 -2.062080  17.0757995            13.095597         41.35950
#chromosome15 25.107664   4.1288069            25.472557         37.87121

rf.bag.prediction<-prediction(rf.bag.pred[,2], data.testing$cancer)
rf.bag.perf<-performance(rf.bag.prediction,"tpr","fpr")
plot(rf.bag.perf,main="ROC Curve for Random Forest (Not caret)",col=2,lwd=2)
abline(a=0,b=1,lwd=2,lty=2,col="gray")
rf.bag.auc <- performance(rf.bag.prediction,"auc")
rf.bag.auc <- unlist(slot(rf.bag.auc, "y.values"))
minauc<-min(round(rf.bag.auc, digits = 2))
maxauc<-max(round(rf.bag.auc, digits = 2))
minauct <- paste(c("min(AUC) = "),minauc,sep="")
maxauct <- paste(c("max(AUC) = "),maxauc,sep="")
minauct
# "min(AUC) = 0.75"
maxauct
# "max(AUC) = 0.75"

#-------------------------
# BUILD THE SVM MODEL - LINEAR KERNAL
#-------------------------
# 36 hours to run
svmFit <- train(cancer ~ ., data = data.training,
                method = "svmLinear",
                trControl = fitControl,
                preProc = c("center", "scale"),
                tuneLength = 8,
                metric = "ROC")
#importance
svmImp <- varImp(svmFit, scale = FALSE)
svmImp
#ROC curve variable importance
#
#  only 20 most important variables shown (out of 29)
#
#             Importance
#seg_size         0.5822
#num_probes       0.5731
#num_genes        0.5721
#most_freq        0.5581
#chromosomeX      0.5077
#chromosome20     0.5055
#chromosome10     0.5052
#chromosome11     0.5050
#chromosome18     0.5047
#chromosome9      0.5045
#chromosome13     0.5045
#chromosome17     0.5045
#chromosome19     0.5037
#chromosome21     0.5035
#chromosome6      0.5015
#chromosome14     0.5012
#chromosome22     0.5007
#chromosome5      0.5004
#chromosome8      0.5003
#chromosome12     0.5002

#-------------------------
# PREDICT USING THE SVM-LINEAR MODEL
#-------------------------
svm.pred<-predict(svmFit,data.testing)
svm.ctable<-confusionMatrix(svm.pred,data.testing$cancer)
svm.ctable
#Confusion Matrix and Statistics
#
#           Reference
#Prediction  cancer no_cancer
#  cancer     12824      3439
#  no_cancer      0         0
#                                          
#               Accuracy : 0.7885          
#                 95% CI : (0.7822, 0.7948)
#    No Information Rate : 0.7885          
#    P-Value [Acc > NIR] : 0.5046          
#                                         
#                  Kappa : 0               
# Mcnemar's Test P-Value : <2e-16          
#                                          
#            Sensitivity : 1.0000          
#            Specificity : 0.0000          
#         Pos Pred Value : 0.7885          
#         Neg Pred Value :    NaN          
#             Prevalence : 0.7885          
#         Detection Rate : 0.7885          
#   Detection Prevalence : 1.0000          
#      Balanced Accuracy : 0.5000          
#                                          
#       'Positive' Class : cancer 
svm.acc<-svm.ctable$overall
svm.acc
#      Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue 
#     0.7885384      0.0000000      0.7821804      0.7947930      0.7885384      0.5045671 
# McnemarPValue 
#     0.0000000 
svmROC<-roc(data.testing$cancer, as.numeric(svm.pred),levels=c("cancer","no_cancer"))
plot(svmROC,main="SVM: 5-fold cross validation")
auc(svmROC)
#Area under the curve: 0.5

#-------------------------
# BUILD THE SVM MODEL - POLYNOMIAL KERNAL
# FULL DATASET
#-------------------------
#2 days to run
svm.model<-train(cancer ~., data=train, method="svmPoly",
                 trControl=trainControl(method='cv',number=5))
svm.model
#Support Vector Machines with Polynomial Kernel 
#
#48788 samples
#   29 predictor
#    2 classes: '1', '-1' 
#
#No pre-processing
#Resampling: Cross-Validated (5 fold) 
#Summary of sample sizes: 39030, 39030, 39031, 39030, 39031 
#Resampling results across tuning parameters:
#
#  degree  scale  C     Accuracy   Kappa      Accuracy SD   Kappa SD   
#  1       0.001  0.25  0.8015906  0.1159665  0.0010104880  0.005505656
#  1       0.001  0.50  0.8021030  0.1208767  0.0009587487  0.005631761
#  1       0.001  1.00  0.8027589  0.1249766  0.0011009689  0.006658491
#  1       0.010  0.25  0.8046036  0.1346137  0.0011371616  0.005492896
#  1       0.010  0.50  0.8054645  0.1411044  0.0008463599  0.005694661
#  1       0.010  1.00  0.8048086  0.1447049  0.0010952670  0.005755187
#  1       0.100  0.25  0.8058539  0.1530854  0.0011949538  0.006246264
#  1       0.100  0.50  0.8060794  0.1559263  0.0009935657  0.006322529
#  1       0.100  1.00  0.8066123  0.1608932  0.0011636536  0.008414941
#  2       0.001  0.25  0.8016521  0.1145189  0.0007314845  0.005018456
#  2       0.001  0.50  0.8023285  0.1188878  0.0009552285  0.004763712
#  2       0.001  1.00  0.8030869  0.1239598  0.0009720548  0.006061570
#  2       0.010  0.25  0.8057105  0.1433472  0.0014170599  0.008045045
#  2       0.010  0.50  0.8076781  0.1570998  0.0012875701  0.007277413
#  2       0.010  1.00  0.8088260  0.1697568  0.0010766577  0.006080419
#  2       0.100  0.25  0.8167377  0.2085598  0.0009967142  0.007275577
#  2       0.100  0.50  0.8185004  0.2201463  0.0006382099  0.005129078
#  2       0.100  1.00  0.8196278  0.2286184  0.0007204266  0.005454871
#  3       0.001  0.25  0.8017751  0.1140491  0.0009079457  0.004535105
#  3       0.001  0.50  0.8022875  0.1185177  0.0009171803  0.004934233
#  3       0.001  1.00  0.8031893  0.1246417  0.0009215288  0.006380319
#  3       0.010  0.25  0.8096253  0.1651637  0.0010410640  0.005766475
#  3       0.010  0.50  0.8113881  0.1775525  0.0009277293  0.007285462
#  3       0.010  1.00  0.8133148  0.1891948  0.0013742494  0.009215494
#  3       0.100  0.25  0.8221899  0.2474767  0.0012894252  0.008947012
#  3       0.100  0.50  0.8233172  0.2551769  0.0011425849  0.008641839
#  3       0.100  1.00  0.8245265  0.2643359  0.0006206245  0.006489139
#
#Accuracy was used to select the optimal model using  the largest value.
#The final values used for the model were degree = 3, scale = 0.1 and C = 1.

#-------------------------
# PREDICT USING THE SVM-POLYNOMIAL MODEL
#-------------------------
svm.poly.pred<-predict(svm.model$finalModel,test[,-30])
svm.poly_ctable<-confusionMatrix(svm.poly.pred,test$cancer)
svm.poly_ctable
svm.poly.acc<-svm.poly_ctable$overall
svm.poly.acc
#Confusion Matrix and Statistics
#
#          Reference
#Prediction     1    -1
#        1  12828  2681
#        -1    69   685
#                                          
#               Accuracy : 0.8309          
#                 95% CI : (0.8251, 0.8366)
#    No Information Rate : 0.793           
#    P-Value [Acc > NIR] : < 2.2e-16       
#                                          
#                  Kappa : 0.2778          
# Mcnemar's Test P-Value : < 2.2e-16       
#                                          
#            Sensitivity : 0.9946          
#            Specificity : 0.2035          
#         Pos Pred Value : 0.8271          
#         Neg Pred Value : 0.9085          
#             Prevalence : 0.7930          
#         Detection Rate : 0.7888          
#   Detection Prevalence : 0.9536          
#      Balanced Accuracy : 0.5991          
#                                          
#       'Positive' Class : 1    

#-------------------------
# BUILD THE SVM MODEL - RADIAL KERNAL
#-------------------------
# 10 days to run
fitControl <- trainControl(method = "repeatedcv",
                           number = 5,
                           repeats = 5,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)
svmFit.rad <- train(cancer ~ ., data = train,
                    method = "svmRadial",
                    trControl = fitControl,
                    #                preProc = c("center", "scale"),
                    tuneLength = 8,
                    metric = "ROC")
RocImp2.svm.rad <- varImp(svmFit, scale = FALSE)
RocImp2.svm.rad
#ROC curve variable importance
#
#  only 20 most important variables shown (out of 29)
#
#             Importance
#seg_size         0.5822
#num_probes       0.5731
#num_genes        0.5721
#most_freq        0.5581
#chromosomeX      0.5077
#chromosome20     0.5055
#chromosome10     0.5052
#chromosome11     0.5050
#chromosome18     0.5047
#chromosome9      0.5045
#chromosome13     0.5045
#chromosome17     0.5045
#chromosome19     0.5037
#chromosome21     0.5035
#chromosome6      0.5015
#chromosome14     0.5012
#chromosome22     0.5007
#chromosome5      0.5004
#chromosome8      0.5003
#chromosome12     0.5002

#-------------------------
# PREDICT USING THE SVM-RADIAL MODEL
#-------------------------
svmRad.pred<-predict(svmFit.rad,data.testing)
svmRad.ctable<-confusionMatrix(svmRad.pred,data.testing$cancer)
svmRad.ctable
#Confusion Matrix and Statistics
#
#           Reference
#Prediction  cancer no_cancer
#  cancer     12698      2730
#  no_cancer    126       709
#                                          
#               Accuracy : 0.8244          
#                 95% CI : (0.8185, 0.8302)
#    No Information Rate : 0.7885          
#    P-Value [Acc > NIR] : < 2.2e-16       
#                                          
#                  Kappa : 0.2716          
# Mcnemar's Test P-Value : < 2.2e-16       
#                                          
#            Sensitivity : 0.9902          
#            Specificity : 0.2062          
#         Pos Pred Value : 0.8230          
#         Neg Pred Value : 0.8491          
#             Prevalence : 0.7885          
#         Detection Rate : 0.7808          
#   Detection Prevalence : 0.9487          
#      Balanced Accuracy : 0.5982          
#                                          
#       'Positive' Class : cancer          
#     
svmROC.rad<-roc(data.testing$cancer, 
                as.numeric(svmRad.pred),levels=c("cancer","no_cancer"))
plot(svmFit.rad,metric="ROC",main="SVM-Radial")
plot(svmROC.rad,main="SVM-Radial: 5-fold cross validation")
auc(svmROC.rad)
#Area under the curve: 0.5982

