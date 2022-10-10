# 10.10.2022
# risk prediction on RNA-seq NB samples using microarray samples used for training NB classifier in Oberthuer et al.
# https://aacrjournals.org/clincancerres/article/21/8/1904/79009/

setwd("/omics/groups/OE0436/data/kurilov/koeln_agilent_custom/2022_rna.seq/")

## 1. RNA-seq data
library(readxl)
library(data.table)
sample_list <- read_excel("NB_evolution_RNAseq.xlsx")
file_list <- paste0(sample_list$Location, sample_list$RNAID_OTP, "/tumor1/paired/merged-alignment/.merging_0/featureCounts/tumor1_", sample_list$RNAID_OTP, ".fpkm_tpm.featureCounts.tsv")
file_list[79] <- gsub("merging_0", "merging_1", file_list[79])

f1 <- fread(file_list[1])
# tpm
tpm_nb <- sapply(file_list, function(x){
  file <- fread(x)
  file$TPM
})

colnames(tpm_nb) <- sample_list$RNAID_OTP
rownames(tpm_nb) <- f1$name
save(tpm_nb, file="tpm_nb.RData")


## 2. Microarray data
library(dplyr)
load("/omics/groups/OE0436/data/kurilov/koeln_agilent_custom/ClinCanRes/rawdata.training.RData")
G.mat <- backgroundCorrect(rawdat.training, method="normexp")
G.mat.quantnorm= normalizeQuantiles(G.mat$E)

ann <- read.delim("GPL16876-29646.txt", skip=25, header=T)
fs_genes <- read.csv("/omics/groups/OE0436/data/kurilov/koeln_agilent_custom/ClinCanRes/features.csv", header = F)
fs_genes <- fs_genes[,1]
fs_genes <- paste0("UKv4_", fs_genes)

fs_num <- match(fs_genes, ann$ProbeName)

exp_microarray <- data.frame(ann$GeneSymbol[fs_num], G.mat.quantnorm[fs_num,])
exp_microarray[,1] <- gsub("C6orf59", "AGPAT4", exp_microarray[,1])
exp_microarray[,1] <- gsub("LOC283454", "HRK", exp_microarray[,1])
exp_microarray2 <- exp_microarray %>% group_by(ann.GeneSymbol.fs_num.) %>% summarise_all(mean)
exp_microarray3 <- data.frame(exp_microarray2[-1,-1])
rownames(exp_microarray3) <- exp_microarray2$ann.GeneSymbol.fs_num.[-1]

exp_training <- t(log2(exp_microarray3))

gene_names_df <- data.frame(old=c("ARPP-21","C19orf36", "C1orf83", "C6orf204",
                                  "C9orf127",  "CRSP2", "FAM36A",   "FAM91A2",  "FLJ22536",
                                  "FLJ22655", "FLJ34503",  "L3MBTL", "LOC199800", "LOC283177",
                                  "LOC286467", "LOC388242", "LOC641784", "RNF12",    
                                  "SKIP",      "SNIP",  "TRA16",  "UGCGL1"),
                            new=c("ARPP21", "IZUMO4", "TCEANC2", "CEP85L",
                                  "TMEM8B", "MED14", "COX20", "LINC00869", "CASC15",
                                  "RERGL", "LINC02880", "L3MBTL1", "ADM5", "B3GAT1",
                                  "FIRRE", "LOC388242", "LOC641784", "RLIM",
                                  "SPHKAP", "SRCIN1", "NR2C2AP", "UGGT1"))
exclude_genes <- setdiff(gene_names_df$new, colnames(exp_test))
sn <- match(gene_names_df$old, colnames(exp_training))
colnames(exp_training)[sn] <- gene_names_df$new
exp_training <- exp_training[,-which(colnames(exp_training) %in% exclude_genes)]

# transforming RNA-seq data
exp_test <- t(log2(tpm_nb+1))
exp_test <- exp_test[,colnames(exp_training)]

# getting class information
load("../ClinCanRes/exprs.quant.training.RData")


## 3. PCA plots
library(ggplot2)
library(ggfortify)
pca_res <- prcomp(exp_training, scale. = TRUE)
autoplot(pca_res, data = exprs.quant.training@phenoData@data, colour = 'class.training')

class_var <- data.frame(class=c(as.character(exprs.quant.training@phenoData@data$class.training), rep(NA, nrow(exp_test))))
pca_res2 <- prcomp(rbind(exp_training, exp_test), scale. = TRUE)
autoplot(pca_res2, data=class_var, colour = 'class')


## 4. Fitting model and generating predictions
library(e1071)
exp_training_s <- scale(exp_training)
exp_test_s <- scale(exp_test)
pca_res3 <- prcomp(rbind(exp_training_s, exp_test_s), scale. = TRUE)
autoplot(pca_res3, data=class_var, colour = 'class')

svm.model <- svm(x=exp_training_s , y=exprs.quant.training@phenoData@data$class.training, cost = 0.015625, gamma = 1, kernel="linear", probability=T, cross=3)
svm.pred <- predict(svm.model, newdata=exp_test_s, probability = T)
df_svm.pred <- data.frame(as.character(svm.pred), attr(svm.pred, "probabilities"))


## 5. Visualizing predictions on PCA plots
class_var2 <- data.frame(class=c(as.character(exprs.quant.training@phenoData@data$class.training), as.character(svm.pred)), batch=rep(c("training", "test"), c(nrow(exp_training_s), nrow(exp_test_s))))
autoplot(pca_res3, data=class_var2, colour = 'class', shape='batch')
autoplot(pca_res2, data=class_var2, colour = 'class')


## 6. Saving predictions
write.csv(cbind(as.character(svm.pred), attr(svm.pred, "probabilities")), file="predictions.csv")
