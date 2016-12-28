# machine-learning-lung-cancer
Project for DSCS6030 (NEU): Comparing Support Vector Machine and Random Forest for Predicting Lung Cancer using Single Nucleotide Polymorphisms

ABSTRACT

Lung cancer risk prediction is particularly difficult in that environmental factors (e.g., smoking) play a huge role in whether or not lung cancer develops. In addition, it is believed that lung cancer is not a disease based on a single gene but is probably associated with the joint effects of many genes. (Yoo, 2012) Both the random forest (RF) and the support vector machine (SVM; Hasibuan, 2014; Fallon, 2013) methods have been shown to be particularly good at identifying multiple predominant features in biological datasets. RF is a classification method that has gained recent popularity for predicting lung cancer risk based on single nucleotide polymorphisms (SNPs). While others have used SVMs to classify SNPs in many different organisms, very few have been applied to classifying SNPs associated with lung cancer. (Hasibuan, 2014; Fallon, 2013) Here, we train two classifiers (SVM and RF) using the class labels of “lung cancer” and “no lung cancer” to identify those features most associated with lung cancer. Our best RF model and SVM models resulted in prediction accuracies of 82.9% and 83.1%, respectively. Both models also determined segment size as the feature with the highest importance. While the RF model had a slightly lower accuracy than the SVM model, creating the RF was much faster and more efficient.

METHODS

~Data Set~

The data were obtained from The Cancer Genome Atlas (https://tcga-data.nci.nih.gov/tcga/) a National Institutes of Health data portal that stores genomic data from tumorous and non-tumorous tissue. It contains a set of individuals where each individual is a collection of SNP data from lung adenocarcinoma tissue along with their corresponding normal tissue samples. There are 110 individuals with each containing a few hundred samples (several per chromosome). There are a total of 220 files: a tumorous sample and normal sample for each individual. The data fields include Sample, Chromosome, Start, End, Num_Probes, and Segment_Mean.

~Preprocessing~

The data were pre-processed in order to create a list of features for training and testing the classifiers. Those features are CHROMOSOME NUMBERS, NUMBER OF PROBES, COPY NUMBER (segment mean), NUMBER OF GENES, VALUE OF THE MOST FREQUENT GENE, and LENGTH OF SEGMENT (segment size). The CHROMOSOME NUMBER is just the number of the chromosome where the SNP is located. The NUMBER OF PROBES is the integer number of probes used in the sequencing experiment. The COPY NUMBER is the log2 ratio of the tumor intensity to the normal intensity (i.e., the segment_mean). The copy number refers to the number of copies of one or more sections of the DNA. As gene duplication occurs sometimes extra copies of a DNA segment are made altering the structure of the DNA. This copy number variation (CNV) results in the cell having an abnormal or, for certain genes, a normal the number of copies of genes. CNV has been linked to cancer progression. CNVs are measured based on the signals of the tissue intensities that are captured during the experiment.  This measurement is the segment mean and represents the log2 ratios of those intensities (Olshen 2004). The values for this feature will be converted to an absolute copy number by the transformation (2^seg_mean)*2. The NUMBER OF GENES is the count of the number genes found in the given SNP segment. The VALUE OF THE MOST FREQUENT GENE is the highest count of the gene most frequently found in that SNP segment. The LENGTH OF SEGMENT is the size of the segment (End position minus Start position).

To calculate the NUMBER OF GENES and VALUE OF MOST FREQUENT GENE, 
parsePos.R parsed each file for chromosome number, and start and end position. These three variables were used to create input files to the UCSC genome browser (http://genome.ucsc.edu/cgi-bin/hgTables). The UCSC genome browser provided output files containing chromosome, gene name, and start and end position of that gene. Using prepDataAddGenes.R, each of these UCSC output files were parsed for gene names, which were appended to the data set to be used for classification. In addition, for each tumorous sample a cancer label of “1” was added, and for each normal sample a label of “-1” was added. Using mergFiles.R, all the preprocessed data files are merged into one .csv file consisting of 65,051 samples with the six features listed above and a label of cancer or no cancer. 

~Classification Models~

Prior to classification, the variables num_probe, num_genes, most_freq, and seg_size were normalized. Training and testing datasets were created based on a 75%-25% split. 

A.  RF model with 5-fold cross validation using the R package caret.

B.  RF without cross validation using the R package randomForest.

C.  SVM/linear kernel with 5-fold cross validation using the R package caret.

D.  SVM/polynomial kernel with 5-fold cross validation using the R package caret.

E.  SVM/radial kernel with 5-fold cross validation using the R package caret.

Two RF models (A and B) were built on the training data set one using 5-fold cross validation and one with out. Model A was trained using the train() function in the R package caret. The trControl parameter was set to “repeatedcv” with 5 repeats. Model B was trained using the randomForest() function in the R package randomForest. Prediction was made on the testing data set using the predict() function from each package. Variable importance was produced using the varImp() function. Calculations of ROC and AUC were done using roc() and auc() functions from the pROC package in R. 

Model B was trained using the randomForest() function from the randomForest R package and predict() was used for the predictions on the testing data set. Variable importance was produced using the importance() function. Calculations of ROC and AUC were done using performance() function from the ROCR package in R.

Three SVM models (C, D, and E) were built on the training data set using 5‐fold cross validation using the train() function in the R package caret. Model C was trained using the method “svmLinear”, which uses a linear kernel. Variable importance was produced using the varImp() function. Prediction was made on the testing data set using the predict() function from the caret package. Calculations of ROC and AUC were done using roc() and auc() functions from the pROC package in R.

Model D was trained using the method “svmPoly”, which uses a polynomial kernel. Prediction was made on the testing data set using the predict() function from the caret package. Variable importance, ROC and AUC were not calculated for model D.

Model E was trained using the method “svmRadial”, which uses a radial kernel. Prediction was made on the testing data set using the predict() function from the caret package. Variable importance was produced using the varImp() function. Calculations of ROC and AUC were done using roc() and auc() functions from the pROC package in R. 

DISCUSSION

~Random Forest~

The two models built using RF had accuracies of 0.8115 and 0.8296 (with 5‐fold cross validation and without, respectively). The higher accuracy of the model not using 5‐fold cross validation is interesting and implies that the model using 5‐fold cross validation may be over-fitted. The RF method builds multiple trees on portions of the data with all the trees averaged at the end (as discussed in "svm_vs_RF_lung_cancer.pdf"). The oob error estimate takes care of this internally. (Breiman, 2001) In addition, the model without 5‐fold cross validation also had better ROC curve and higher AUC (0.75 versus 0.6234).

Both of these models listed chromosome 6 and chromosome 7 as variables of importance after seg_size, num_probes, and num_genes. Since all models created here (both RF and SVM) listed seg_size, num_probes, and num_genes as variables of importance, using chromosomes as features for prediction may not be ideal.

~Support Vector Machine~

The three models (C, D, and E) built using SVM had accuracies of 0.7885, 0.8309, and 0.8244 (with linear, polynomial, and radial kernels, respectively). Model C classified all observations as “cancer”, and by doing so resulted in a fairly high accuracy.
However, this model produced an ROC and an AUC (0.5) that are essentially identical to what would be produced through random guessing.

While the SVM‐polynomial model had the highest accuracy, it took appreciably more time (days versus hours) to create than the SVM‐linear or either of the RF models.
While the SVM‐polynomial model had the highest accuracy, it is too time-consuming. In addition, due to the long computing times, the importance variables, ROC and AUC metrics were not produced. Without these, it is difficult to determine whether this higher accuracy is meaningful.

Model E took, by far, the longest time to train at 10 days. This model and Model B (RF without cross validation) had very comparable results, e.g., similar sensitivities, specificities. In addition, each of their accuracies were within 95% CI of one another. Due to the similarities of the results of these two models and the large discrepancy in their computation time, Model B is the preferred model here. As a follow-up, it would be interesting to see whether Model E’s training time could be reduced enough through distributed computation scheme that would make it the preferred model over any RF model.

For additional details, please see the pdf "svm_vs_RF_lung_cancer.pdf".

REFERENCES

Boser, B. E., Isabelle M. Guyon , Vladimir N. Vapnik. A Training Algorithm for Optimal Margin Classifiers. Proceedings of the 5th Annual ACM Workshop on Computational Learning Theory. 1992.

Breiman, Leo. "Random Forests". Machine Learning 45 (1): 5–32. 2001.

Cortes, C. and Vapnik, V. Support-­‐Vector Networks. Machine Learning, 20, 273-­‐297. 1995.

Fallon, B. O., Whitney Wooderchak-­‐Donahue, and David K. Crockett, A support vector machine for identification of single-­‐nucleotide polymorphisms from next-­‐ generation sequencing data, Bioinformatics, Vol. 29 no. 11 2013, pages 1361–1366

Fawcett, T. An introduction to ROC analysis. Pattern Recognition Letters 27 (2006) 861–874

Hasibuan, L.S. Kusuma, W.A. ; Suwamo, W.B. Identification of single nucleotide polymorphism using support vector machine on imbalanced data. Advanced Computer Science and Information Systems (ICACSIS), 2014 International Conference. 375 – 379.

Rong Li, Jing Wu, Liwen Xiong, and Baohui Han. Lung cancer and benign lung diseases in patients with serious vitamin D deficiency in eastern China. Thoracic Cancer Volume 3, Issue 4, pages 303–306, November 2012.

Teknomo, Kardi. (2009) Tutorial on Decision Tree. http://people.revoledu.com/kardi/tutorial/DecisionTree/

Yoo, W., et al., A Comparison of Logistic Regression, Logic Regression, Classification Tree, and Random Forests to Identify Effective Gene-­‐Gene and Gene-­‐Environ.
International Journal of Applied Science and Technology Vol. 2 No. 7; August 2012.


