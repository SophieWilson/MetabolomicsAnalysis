if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
biocLite()
BiocManager::install("impute")
BiocManager::install("ProteoMM")
library(ProteoMM)
library(impute)
setwd('C:/Users/Mischa/Documents/Uni Masters/Module 4')
getwd()
Impute_Normalise <- read.table('Galaxy6-[Normalisation_on_data_3_and_data_5__Peak_intensity_matrix_(Imputation_then_normalisation)].tsv', header = TRUE)
Normalise_impute <- read.table('Galaxy7-[Missing_values_imputation_on_data_3_and_data_4__Peak_intensity_matrix_(normalisation_then_imputation)].tsv', header = TRUE)
ID_df <- read.table('Galaxy3-[Galaxy16-[Blank_Filter_on_data_10__Sample_Metadata_(updated)].tsv].tabular', header = TRUE)
ID_df <- ID_df[-c(1,7,13,19,25),]


Normalise_impute <- Normalise_impute[-c(1,7,13,19,25),]
Impute_Normalise <- Impute_Normalise[-c(1,7,13,19,25),]
group1 <- cbind(ID_df, Impute_Normalise)
group1 <- group1[,-c(2,3,4,6)]
group2 <- cbind(ID_df, Normalise_impute)
write.table(group1, "group1.txt")

levels(group1$classLabel)
group1 <- with(group1, group1[order(classLabel),])
group2 <- with(group2, group2[order(classLabel),])

# pca for both groups
pca_group1 <- Impute_Normalise[,-1]
pca_resgrp1 <- prcomp(pca_group1, center = TRUE, scale = TRUE)
pc_df1 <- data.frame(pca_resgrp1$x, Classlabel = group1$classLabel, 
                     name = group1$filename)
ggplot(pc_df1, aes(x = PC1, y = PC2, shape = Classlabel, col = name)) +
  geom_point(size = 6)
             
pca_group2 <- Normalise_impute[,-1]
pca_resgrp2 <- prcomp(pca_group2, center = TRUE, scale = TRUE)
pc_df2 <- data.frame(pca_resgrp2$x, Classlabel = group2$classLabel, 
                     name = group2$filename)
ggplot(pc_df2, aes(x = PC1, y = PC2, shape = Classlabel)) +
  geom_point(size = 6)             
             
summary(pca_resgrp1)
summary(pca_resgrp2)

# making some covariance matrices 
cov_1 <- cov(pca_group1)
cov_2 <- cov(pca_group2)

# total variance of each dataframe, the sum of eigenvalues
sum(diag(cov_1))
sum(diag(cov_2))

s.eigen1 <- eigen(cov_1)
s.eigen2 <- eigen(cov_2)
sum(s.eigen1$values[1:2])
sum(s.eigen2$values[1:2])
for (s in 1:2){
  print(s/sum(s.eigen$values))
}



plot(s.eigen2$values[1:6], xlab = 'Eigenvalue Number', ylab = 'Eigenvalue Size', main = 'Scree Graph')
lines(s.eigen2$values[1:6])

res.aov <- aov(weight ~ group, data = my_data)

# removing quality control samples
ID_df2 <- ID_df[-c(1,7,13,19,25),]

ID_df2$secondclass <- ID_df2$classLabel
# making the labels in the correct format
metaCols = 5:6 # column indices for metadata such as protein IDs and sequences
m_prot.info = make_meta(ID_df2, metaCols)
prot.info <- cbind(ID_df2$classLabel, ID_df2$classLabel)
# making the groups numerical 
groups <- as.factor(ID_df2$classLabel)
# integer with some values
intscols <- 2-30
m_logInts = make_intencities(Normalise_impute[,-1])  
m_logInts = convert_log2(m_logInts) 
plot = eig_norm1(m = m_logInts,treatment=groups,prot.info=m_prot.info)
print(groups)


####################################################################
set.seed(1234)
df_1 <- cbind(ID_df, Impute_Normalise[,-1])
df_1 <- df_1[-c(1,7,13,19,25),]
df_1 <- t(df_1)
df_1 <- df_1[-c(1,2,3,4),]
colnames(df_1) <- df_1[1,]
df_1 <- df_1[-1,]
metaCols = 1:5
intsCols <- 6:12
m_logInts = data.frame(make_intencities(df_1)) # will reuse the name
dim(m_logInts)
m_logInts <- m_logInts[ , order(names(m_logInts))]
m_prot.info = data.frame(cbind(rownames(df_1), rownames(df_1)))
m_prot.info <- make_meta(m_prot.info)
dim(m_prot.info)
# 3 samples for CG and 3 for mCG
groups <- as.factor(c('cow', 'cow', 'cow', 'cow', 'cow', 'cow', 'cow', 'cow', 'cow', 'cow', 'sheep', 'sheep', 'sheep', 'sheep', 'sheep', 'sheep','sheep', 'sheep', 'sheep', 'sheep'))
length(groups)
mm_m_ints_eig1 = eig_norm1(m=m_logInts,treatment=groups,prot.info=m_prot.info)


###############################################################
# From source code, this works but have to adjust dataset. 
data("mm_peptides")
head(mm_peptides)
# different from parameter names as R uses outer name spaces
# if variable is undefined
intsCols = 8:13
metaCols = 1:7 # reusing this variable
m_works = make_intencities(mm_peptides, intsCols) # will reuse the name
dim(m_works)
info_works = make_meta(mm_peptides, metaCols)
dim(info_works)
m_works = convert_log2(m_works)
# 3 samples for CG and 3 for mCG
grps = as.factor(c('CG','CG','CG', 'mCG','mCG','mCG'))
length(grps)
# ATTENTION: SET RANDOM NUMBER GENERATOR SEED FOR REPRODUCIBILITY !!
set.seed(123) # Bias trends are determined via a permutaion, results may
# vary slightly if a different seed is used, such as when set.seed()
# function is not used
mm_m_ints_eig1 = eig_norm1(m=m_works,treatment=grps,prot.info=info_works)
mm_m_ints_eig1$h.c # check the number of bias trends detected
mm_m_ints_norm = eig_norm2(rv=mm_m_ints_eig1)


unprocessed_dat <- read.csv('batch7_SFPM_unprocessed (1).csv')
getwd()
setwd('C:/Users/Mischa/Documents/Uni Masters/Module 4')

install.packages('missForest')
install.packages('DMwR')
library(missForest)
library(DMwR)

#'' this uses normalised data to predict missing values 
unprocessed <- read.table('Galaxy20-[Normalisation_on_data_16_and_data_18__Peak_intensity_matrix] (3).tsv', header = TRUE)

with_header <- unprocessed
unprocessed_rm <- unprocessed[,-c(1)]
unprocessed_rm[unprocessed_rm == 0] <- NA
processed <- missForest(unprocessed_rm)
processeddf <- processed$ximp




## testing the accuracy of the RF imputation ### 

write.csv(processeddf, 'C:/Users/Mischa/Documents/Uni Masters/Module 4/normalised_rfimputed.csv' )
## first import all the data i can and will merge the ones with the same columns
batch_4 <- read.table('galaxy_normalisedfulldata.tsv', header = TRUE)
batch_5 <- read.table('batch5.tsv', header = TRUE)
batch_6 <- read.table('batch6.tsv', header = TRUE)
batch_7 <- read.table('batch7.tsv', header = TRUE)
batch_9 <- read.table('batch9.tsv', header = TRUE)

# getting as much data as possible
full <- rbind(batch_6,batch_5)
full[full == 0] <- NA
# removing any columns with na (every row has na so cant do it that way)
no_natest <- full[ , colSums(is.na(full)) == 0]
no_natest <- no_natest[,-1]
# testing the random forest imputation

# adding 15% of missing values
added_na <- as.data.frame(lapply(no_natest, function(cc) cc[ sample(c(TRUE, NA), prob = c(0.85, 0.15), size = length(cc), replace = TRUE) ]))
## imputing the forset, and testing against true values
full_impute <- missForest(added_na, xtrue = no_natest)
full_impute$OOBerror

####### someone elses graph ###########
train_test_X <- added_na
imp_test_X <- missForest(train_test_X)$ximp[1:nrow(test_X), ]

impute.knn(train_test_X, k=10)


# calculating error
sum_true <- colSums(no_natest)
sum_predict <- colSums(full_impute$ximp)
difference <- 0
error <- 0
for (x in 1:length(sum_true)) {
  difference[x] <- (sum_true[x] - sum_predict[x])
  error[x] <- difference[x] / sum_true[x]
  
}
sum(error)
mean(difference)

install.packages('HDclassif')
library(HDclassif)
hddc <- hddc(processeddf, K = 2)

### using KNN imputation
# knnOutput <- knnImputation(unprocessed_rm[, !names(unprocessed_rm) %in% "medv"], k=1)  # perform knn imputation.
# anyNA(knnOutput)

library(class)
pr <- knn(unprocessed_rm,k=2)




########### looking at different features of cow and sheep #########
library(ggplot2)


labels <- read.csv('batch7_SFPM_fully_processed (1).csv', header = TRUE)
processeddf$labels <- labels$Label
pca_predata <- processeddf[processeddf$labels != "QC", ]
nolabel <- pca_predata[,-677]
pca_res1 <- prcomp(nolabel, center = TRUE, scale = TRUE)
pc_df1 <- data.frame(pca_res1$x, class = pca_predata$labels)
ggplot(pc_df1, aes(x = PC1, y = PC2, col = class)) +
  geom_point(size = 6) + 
  ggtitle(~underline('PCA plot for patient gene expression')) +
  theme(plot.title = element_text(hjust=0.5)) 



write.csv(pca_predata, 'C:/Users/Mischa/Documents/Uni Masters/Module 4/pcadf.csv' )
########### using someone elses code #########
install.packages("corrplot")
install.packages('factoextra')
library("corrplot")
library("factoextra")

var <- get_pca_var(pca_res1)
corrplot(var$cos2, is.corr=FALSE)
fviz_cos2(pca_res1, choice = "var", axes = 1:2)
fviz_contrib(pca_res1, choice = "ind", axes = 1:2)

head(var$contrib, 4)
#Call the function. Use only the 2 PCs.
myplot(x_new[:,0:2],np.transpose(pca.components_[0:2, :]))
plt.show()

fviz_pca_biplot(pca_res1, select.ind = list(contrib = 5), 
                select.var = list(contrib = 5),
                ggtheme = theme_minimal())

for (x in 1:length(var$cos2)){
  
}


library(class)
##run knn function


##create confusion matrix
tab <- table(pr,iris_test_category)

##this function divides the correct predictions by total number of predictions that tell us how accurate teh model is.

accuracy <- function(x){sum(diag(x)/(sum(rowSums(x)))) * 100}
accuracy(tab)
