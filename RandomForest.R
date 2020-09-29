##############################
#### 1) Required packages ####
##############################

if(!require("openxlsx")) (install.packages("openxlsx")); library(openxlsx)
if(!require("readr")) (install.packages("readr")); library(readr)
if(!require("caret")) (install.packages("caret")); library(caret)
if(!require("kableExtra")) (install.packages("kableExtra")); library(kableExtra)
if(!require("ggplot2")) (install.packages("ggplot2")); library(ggplot2)
if(!require("ROCR")) (install.packages("ROCR")); library(ROCR)
if(!require("pROC")) (install.packages("pROC")); library(pROC)
if(!require("gplots")) (install.packages("gplots")); library(gplots)
if(!require("RColorBrewer")) (install.packages("RColorBrewer")); library(RColorBrewer)
if(!require("nnet")) (install.packages("nnet")); library(nnet)
if(!require("randomForest")) (install.packages("randomForest")); library(randomForest)
if(!require("caret")) (install.packages("caret")); library(caret)
if(!require("VIM")) (install.packages("VIM")); library(VIM)
if(!require("laeken")) (install.packages("laeken")); library(laeken)
if(!require("ggbiplot")) (install.packages("ggbiplot")); library(ggbiplot)
if(!require("dplyr")) (install.packages("dplyr")); library(dplyr)
if(!require("reshape")) (install.packages("reshape")); library(reshape)
# For the decision tree:
library(devtools)
install_github("araastat/reprtree")
library(reprtree)
if(!require("gridExtra")) (install.packages("gridExtra")); library(gridExtra)
if(!require("grid")) (install.packages("grid")); library(grid)

########################
#### 2) Source data ####
########################

source("read_ddbb_RF.R")
summary(df)

################################ 
#### 3) Imput NA's with kNN ####
################################ 

newdata <- VIM::kNN(df, variable = c("blood_white_cell_count_t0", "CRP_t0", "Hemoglobin_t0", "number_relapses", "Disease_duration"), k = 5, numFun = weightedMean, weightDist=TRUE)
df <- newdata[,1:24]
#summary(df)

##############################################
#### 4) Split data in train and test sets ####
##############################################

df$Disease <- factor(df$Disease, levels = c("UC", "CD"))
contrasts(df$Disease)

set.seed(1234)
indexes <- createDataPartition(y = df$Disease, p = 2/3, list = FALSE)
train <- df[indexes,]
test <- df[-indexes,]

table(train$Disease)
table(test$Disease)

###################################################################################
#### 5) Function of 72 random forest combinations for different mtry and ntree ####
###################################################################################

rf_models <- function(explicative.vars, ntrees, mtrys){
  vars <- paste(explicative.vars, collapse="+") 
  result <- c()
  ids <- c()
  for(i in ntrees){
    for(j in mtrys){

        #set.seed(123456)
        classifier1 <- randomForest(as.formula(paste("Disease ~ ", vars, sep = "")),
                                    data=train, ntree=i, mtrys=j, importance=TRUE)
        
        badclassified_train <- train[classifier1$predicted!=classifier1$y,]
        badclassified_train_ids <- badclassified_train[,"ID"]
        badclassified_train_ids <- badclassified_train_ids[!is.na(badclassified_train_ids)]
        
        pred <- predict(classifier1, test)
        badclassified <- test[,"Disease"]!=pred
        badclassified_test_ids <- test[badclassified,"ID"]
        
        ids <- c(ids, badclassified_test_ids, badclassified_train_ids)
        
        conf_mat <- confusionMatrix(pred, as.factor(test$Disease))
        accuracy <- conf_mat$overall[1]
        sensitivity <- conf_mat$byClass[1]
        specificity <- conf_mat$byClass[2]
        r <- c(i, j,  accuracy, sensitivity, specificity)
        result <- rbind(result, r)
      
    }
  }
  colnames(result) <- c("Tree", "Variables", "Accuracy", "Sensitivity", "Specificity")
  return(list(result,ids))
}

##################################
#### 6) Evaluate the function ####
##################################

# Fix different values for ntree and mtry:
  
ntrees <- c(5, 10, 20, 50, 100, 200, 500, 1000)
mtrys <- c(1, 2, 3, 5, 7, 9, 12, 15, 19) 

# Select explicative variables and apply the function:
  
explicative.vars <- c("ITS2_t0", "V4_t0", "Ratio_t0", "Gender", "Age", "Height_m", "Weight_Kg", "BMI", 
                      "Smoking", "Calprotectin_t0", "blood_white_cell_count_t0", "CRP_t0", "Hemoglobin_t0", 
                      "number_relapses", "Disease_duration", "Shannon", "Chao1", "Escherichia", "Clostridium",
                      "Coprobacillus1", "Faecalibacterium")
result <- rf_models(explicative.vars, ntrees, mtrys)

#### 6.1) Barplot of the missclassified individuals
badclassified_ids <- result[[2]]
badclassified <- as.data.frame(table(badclassified_ids))
badclassified$badclassified_ids_corrected <- substr(badclassified$badclassified_ids, 5, 15)
#pdf("Barplot_CDUC_badclassified.pdf") 
ggplot(badclassified, aes(x = reorder(badclassified_ids_corrected, -Freq), y = Freq)) + 
  geom_bar(stat = "identity") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 0.5, size=12)) +
  labs(x="") + 
  theme(axis.text = element_text(size=14, face = "bold")) +
  theme(axis.title.y = element_text(size=16, face = "bold")) 
#dev.off()

#### 6.2) Cluster of the missclassified individuals
newdf <- df
rownames(newdf) <- newdf[,1]
newdf <- subset(newdf, newdf$ID %in% badclassified$badclassified_ids)
newdf$ID <- substr(newdf$ID, 5, 15)
rownames(newdf) <- substr(rownames(newdf), 5, 15)
# Ward Hierarchical Clustering
d <- dist(newdf[,-c(1,2)], method = "euclidean") # distance matrix
fit <- hclust(d, method="ward") 
#pdf("Cluster_CDUC_badclassified.pdf") 
plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit, k=5, border="red")
#dev.off()

########################################################################
#### 7) Mean decrease Gini and mean decrease accuracy in all models ####
########################################################################
result2 <- as.data.frame(result[[1]])
head(result2)

ntree <- as.numeric(paste(result2$Tree)) # all combinations for ntree and mtry from above
mtry <- as.numeric(paste(result2$Variables))

explicative.vars <- c("ITS2_t0", "V4_t0", "Ratio_t0", "Gender", "Age", "Height_m", "Weight_Kg", "BMI", 
                      "Smoking", "Calprotectin_t0", "blood_white_cell_count_t0", "CRP_t0", "Hemoglobin_t0", 
                      "number_relapses", "Disease_duration", "Shannon", "Chao1", "Escherichia", "Clostridium",
                      "Coprobacillus1", "Faecalibacterium")
vars <- paste(explicative.vars, collapse="+") 

data_MDAccuracy <- c()
data_MDGini <- c()

for(i in 1:length(ntree)){
  # set.seed(123456)
  classifier <- randomForest(as.formula(paste("Disease ~ ", vars, sep = "")),
                             data=train, ntree=ntree[i], mtrys=mtry[i], importance=TRUE)
  
  df_importance <- randomForest::importance(classifier)
  
  MeanDecreaseAccuracy <- as.data.frame(df_importance[,3])
  colnames(MeanDecreaseAccuracy) <- paste0("ntree=",ntree[i], ";mtry=",mtry[i])
  
  MeanDecreaseGini <- as.data.frame(df_importance[,4])
  colnames(MeanDecreaseGini) <- paste0("ntree=",ntree[i], ";mtry=",mtry[i])
  # Create data for boxplots
  data_MDAccuracy <- cbind(data_MDAccuracy, MeanDecreaseAccuracy[,1])
  colnames(data_MDAccuracy)[i] <- colnames(MeanDecreaseAccuracy)[1]
  rownames(data_MDAccuracy) <- rownames(MeanDecreaseAccuracy)
  
  data_MDGini <- cbind(data_MDGini, MeanDecreaseGini[,1])
  colnames(data_MDGini)[i] <- colnames(MeanDecreaseGini)[1]
  rownames(data_MDGini) <- rownames(MeanDecreaseGini)
}

newdf <- t(data_MDAccuracy)
newdf <- melt(newdf)
#pdf("MeanDecreaseAccuracy.pdf") 
p1 <- ggplot(newdf, aes(x = X2, y = value)) +
  geom_boxplot() + 
  theme_bw() + 
  scale_x_discrete(labels=c("blood_white_cell_count_t0" = "blood_white_\ncell_count_t0")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size=11)) + 
  theme(axis.text = element_text(size=14, face = "bold")) + 
  labs(x="", y="") + 
  theme(plot.title = element_text(face = "bold", size=14)) + 
  ggtitle("Mean Decrease Accuracy")
print(p1)
#dev.off()

newdf <- t(data_MDGini)
newdf <- melt(newdf)
#pdf("MeanDecreaseGini.pdf") 
p2 <- ggplot(newdf, aes(x = X2, y = value)) +
  geom_boxplot() + 
  theme_bw() + 
  scale_x_discrete(labels=c("blood_white_cell_count_t0" = "blood_white\ncell_count_t0")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size=11)) + 
  theme(axis.text = element_text(size=14, face = "bold")) + 
  labs(x="", y="") + 
  theme(plot.title = element_text(face = "bold", size=14)) + 
  ggtitle("Mean Decrease Gini")
print(p2)
#dev.off()

#grid.arrange(p1, p2, ncol = 2)


########################################
#### 8) Hyperparameter optimization ####
########################################

#### 8.1) Identification of the optimal value of the hyperparameter mtry.

tuning_rf_mtry <- function(df_train, explicative.vars, ntree = 500){
  # Output = out-of-bag clasification error 
  # explicative.vars = predictor variables
  
  require(dplyr)
  vars <- paste(explicative.vars, collapse="+") 
  max_predictores <- length(explicative.vars)
  n_predictores   <- rep(NA, max_predictores)
  oob_err_rate    <- rep(NA, max_predictores)
  for (i in 1:max_predictores) {
    set.seed(123)
    
    model_rf <- randomForest(as.formula(paste("Disease ~ ", vars, sep = "")), 
                             data = df_train, mtry = i, ntree = ntree)
    n_predictores[i] <- i
    oob_err_rate[i] <- tail(model_rf$err.rate[, 1], n = 1)
  }
  results <- data_frame(n_predictores, oob_err_rate)
  return(results)
}

hiperparametro_mtry <-  tuning_rf_mtry(df_train=train, explicative.vars)
hiperparametro_mtry %>% arrange(oob_err_rate) %>% head()

p3 <- ggplot(data = hiperparametro_mtry, aes(x = n_predictores, y = oob_err_rate)) +
  scale_x_continuous(breaks = hiperparametro_mtry$n_predictores) +
  geom_line() +
  geom_point() +
  theme_bw() + 
  geom_point(data = hiperparametro_mtry %>% arrange(oob_err_rate) %>% head(1),
             color = "red") +
  labs(title = "Out-of-bag-error vs mtry",
       x = "mtry") +
  theme(axis.text = element_text(size=8, face = "bold")) + 
  theme(axis.title = element_text(size=11, face = "bold")) + 
  theme(plot.title = element_text(face = "bold", size=12)) 
print(p3)

best_mtry <- which.min(hiperparametro_mtry$oob_err_rate)


#### 8.1.2) Identification of the optimal value of the nodesize hyperparameter

tuning_rf_nodesize <- function(df_train, size = NULL, ntree = 500){
  # Output = out-of-bag clasification error
  # df_train = train set
  # sizes = evaluated sizes
  
  require(dplyr)
  if (is.null(size)){
    size <- seq(from = 1, to = nrow(df), by = 5)
  }
  oob_err_rate <- rep(NA, length(size))
  vars <- paste(explicative.vars, collapse="+") 
  
  for (i in seq_along(size)) {
    set.seed(321)
    
    model_rf <- randomForest(as.formula(paste("Disease ~ ", vars, sep = "")), data = df_train, 
                             mtry = best_mtry, ntree = ntree, nodesize = i)
    oob_err_rate[i] <- tail(model_rf$err.rate[, 1], n = 1)
  }
  results <- data_frame(size, oob_err_rate)
  return(results)
}

hiperparametro_nodesize <-  tuning_rf_nodesize(df = train,
                                               size = c(1:20))
hiperparametro_nodesize %>% arrange(oob_err_rate) %>% head()

p4 <- ggplot(data = hiperparametro_nodesize, aes(x = size, y = oob_err_rate)) +
  scale_x_continuous(breaks = hiperparametro_nodesize$size) +
  geom_line() +
  geom_point() +
  theme_bw() + 
  geom_point(data = hiperparametro_nodesize %>% arrange(oob_err_rate) %>% head(1),
             color = "red") +
  labs(title = "Out-of-bag-error vs nodesize",
       x = "nodesize") +
  theme(axis.text = element_text(size=8, face = "bold")) + 
  theme(axis.title = element_text(size=11, face = "bold")) + 
  theme(plot.title = element_text(face = "bold", size=12)) 
print(p4)

best_nodesize <- which.min(hiperparametro_nodesize$oob_err_rate)

#### 8.1.3) Identification of the optimal value of the hyperparameter ntree

vars <- paste(explicative.vars, collapse="+") 
modelo_randomforest <- randomForest(as.formula(paste("Disease ~ ", vars, sep = "")), data = train, mtry = best_mtry, ntree = 1000, importance = TRUE, nodesize = best_nodesize)
oob_err_rate <- data.frame(oob_err_rate = modelo_randomforest$err.rate[, 1],
                           arboles = seq_along(modelo_randomforest$err.rate[, 1]))
p5 <- ggplot(data = oob_err_rate, aes(x = arboles, y = oob_err_rate )) +
  geom_line() +
  theme_bw() + 
  labs(title = "Out-of-bag-error vs ntree",
       x = "ntree") +
  theme(axis.text = element_text(size=8, face = "bold")) + 
  theme(axis.title = element_text(size=11, face = "bold")) + 
  theme(plot.title = element_text(face = "bold", size=12))
print(p5)

best_ntree <- 500


###########################################
#### 9) Final fit and out-of-bag error ####
###########################################

#### 9.1) Final fit and out-of-bag (OOB) error

explicative.vars <- c("ITS2_t0", "V4_t0", "Ratio_t0", "Gender", "Age", "Height_m", "Weight_Kg", "BMI", "Smoking", "Calprotectin_t0", "blood_white_cell_count_t0", "CRP_t0", "Hemoglobin_t0", "number_relapses", "Disease_duration", "Shannon", "Chao1", "Escherichia", "Clostridium", "Coprobacillus1", "Faecalibacterium")
vars <- paste(explicative.vars, collapse="+") 
allvars_model <- randomForest(as.formula(paste("Disease ~ ", vars, sep = "")), data = train, mtry = best_mtry, 
                              ntree = best_ntree, importance = TRUE, nodesize = best_nodesize, norm.votes = TRUE )
allvars_model


explicative.vars <- c("Ratio_t0", "blood_white_cell_count_t0", "Hemoglobin_t0", "Shannon", "Chao1", "Escherichia", "Clostridium", "Faecalibacterium")
vars <- paste(explicative.vars, collapse="+") 
final_model <- randomForest(as.formula(paste("Disease ~ ", vars, sep = "")), data = train, mtry = best_mtry, 
                            ntree = best_ntree, importance = TRUE, nodesize = best_nodesize, norm.votes = TRUE, proximity=TRUE )
final_model

# Difference in error classification using all the explicative variables or only the most important:
# Last value OOB
final_model_oob <- tail(final_model$err.rate[, 1], 1)
final_model_oob
# Last value OOB
allvars_model_oob <- tail(allvars_model$err.rate[, 1], 1)
allvars_model_oob


#### 9.2) Some plots of the final classifier to help its interpretation:

# Margin: positive margin means correct classification, and vice versa.

pred <- predict(final_model, test)
badclassified <- test[,"Disease"]!=pred
badclassified_id <- test[badclassified,"ID"]
table(badclassified_id)
plot(margin(final_model))
abline(h=0)

# MDSplot

MDSplot(final_model, fac=train$Disease)

# partiaPlot: gives a graphical depiction of the marginal effect of a variable on the class probability 
#             (classification) or response (regression).

library(devtools)
devtools::install_github("zmjones/party", subdir = "pkg")
devtools::install_github("zmjones/edarf", subdir = "pkg")
library(edarf)

# Get variable importance measures
imp_df <- data.frame(importance(final_model, scale = FALSE, type = 1))
# Tidy up and sort the data frame
imp_df <- imp_df %>% 
  mutate(names = rownames(imp_df)) %>% 
  arrange(desc(MeanDecreaseAccuracy))
# Plot mean decreased accuracy
imp_df %>% 
  top_n(8, MeanDecreaseAccuracy) %>% 
  ggplot(aes(x = reorder(names, MeanDecreaseAccuracy),y = MeanDecreaseAccuracy)) +
  geom_col() +
  coord_flip() +
  labs(title = "Variable Importance",
       subtitle = "Random Forests",
       x= "",
       y= "Mean Decrease in Accuracy")

# Save top predictor names as character vector
nm <- as.character(imp_df$names)[1:8]
# Get partial depedence values for top predictors
pd_df <- partial_dependence(fit = final_model,
                            vars = nm,
                            data = train)
#pdf("Partial_dependence_plots.pdf") 
plot_pd(pd_df) +
  theme_bw() + 
  scale_color_manual(name = "Disease",
                     values = c("#5F7FC7", "orange")) + 
  theme(strip.text = element_text(face = "bold", size = 12),
        legend.title = element_text(size=13, face="bold"),
        legend.text = element_text(size=12, face="bold")) +
  labs(x="") 
#dev.off()

#### 9.3) Most influential predictors using the best values for the hyperparameters

if(require("tidyverse")) (install.packages("tidyverse")); library(tidyverse)
if(require("ggpubr")) (install.packages("ggpubr")); library(ggpubr)
if(require("pBrackets")) (install.packages("pBrackets")); library(pBrackets)

importancia_pred <- as.data.frame(importance(allvars_model, scale = TRUE))
importancia_pred <- rownames_to_column(importancia_pred, var = "variable")

q1 <- ggplot(data=importancia_pred, aes(x=reorder(variable, MeanDecreaseAccuracy),
                                        y = MeanDecreaseAccuracy)) +
  labs(x="") +
  geom_col() +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom") + 
  theme(axis.text = element_text(size=9, face="bold")) +
  theme(axis.title = element_text(size=9, face = "bold"))

q2 <- ggplot(data = importancia_pred, aes(x = reorder(variable, MeanDecreaseGini),
                                          y = MeanDecreaseGini)) +
  labs(x = "") +
  geom_col() +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom") + 
  theme(axis.text = element_text(size=9, face="bold")) + 
  theme(axis.title = element_text(size=9,face = "bold")) 

grid.arrange(q1, q2, ncol = 2)


##########################################
#### 10) ROC Curve and decision trees ####
##########################################

#### 10.1) ROC Curve

###############
explicative.vars <- c("blood_white_cell_count_t0", "Hemoglobin_t0")
vars <- paste(explicative.vars, collapse="+") 

classifier1 <- randomForest(as.formula(paste("Disease ~ ", vars, sep = "")), data=train, ntree=best_ntree, mtrys=best_mtry, nodesize = best_nodesize,importance=TRUE)

predictions1 <- as.data.frame(predict(classifier1, test, type = "prob"))
predictions1$predict <- names(predictions1)[1:2][apply(predictions1[,1:2], 1, which.max)]
predictions1$observed <- test$Disease
roc1 <- roc(ifelse(predictions1$observed=="CD", "CD", "UC"), as.numeric(predictions1$CD))
pROC::coords(roc1, x="best", ret=c("specificity", "sensitivity"), input="threshold", best.method="youden") 

###############
###############
explicative.vars <- c("Ratio_t0", "Shannon", "Chao1", "Escherichia", "Clostridium", "Faecalibacterium")
vars <- paste(explicative.vars, collapse="+") 

classifier2 <- randomForest(as.formula(paste("Disease ~ ", vars, sep = "")), data=train, ntree=best_ntree, mtrys=best_mtry, nodesize = best_nodesize,importance=TRUE)

predictions2 <- as.data.frame(predict(classifier2, test, type = "prob"))
predictions2$predict <- names(predictions2)[1:2][apply(predictions2[,1:2], 1, which.max)]
predictions2$observed <- test$Disease
roc2 <- roc(ifelse(predictions2$observed=="CD", "CD", "UC"), as.numeric(predictions2$CD))
pROC::coords(roc2, x="best", ret=c("specificity", "sensitivity"), input="threshold", best.method="youden") 

###############
###############
explicative.vars <- c("Ratio_t0", "blood_white_cell_count_t0", "Hemoglobin_t0", "Shannon", "Chao1", "Escherichia", "Clostridium", "Faecalibacterium")
vars <- paste(explicative.vars, collapse="+") 

classifier3 <- randomForest(as.formula(paste("Disease ~ ", vars, sep = "")), data=train, ntree=best_ntree, mtrys=best_mtry, nodesize = best_nodesize,importance=TRUE)

predictions3 <- as.data.frame(predict(classifier3, test, type = "prob"))
predictions3$predict <- names(predictions3)[1:2][apply(predictions3[,1:2], 1, which.max)]
predictions3$observed <- test$Disease
roc3 <- roc(ifelse(predictions3$observed=="CD", "CD", "UC"), as.numeric(predictions3$CD))
pROC::coords(roc3, x="best", ret=c("specificity", "sensitivity"), input="threshold", best.method="youden") 

##############
#### PLOT ####
##############
q3 <- ggroc(list(Combination = roc3, Clinical = roc1, Microbiota=roc2), aes = c("linetype","color")) +
  geom_line(size=1) +
  geom_abline(intercept = 1, colour="gray") + 
  theme_bw() + 
  scale_color_manual(values=c("darkgreen","darkblue","yellow"),
                     labels=c(paste0("Combination, AUC=",round(roc3$auc,3)), paste0("Clinical, AUC=", round(roc1$auc,3)), paste0("Microbiota, AUC=", round(roc2$auc,3)))) +
  scale_linetype_manual(values = c(1,2,3),
                        labels=c(paste0("Combination, AUC=",round(roc3$auc,3)), paste0("Clinical, AUC=", round(roc1$auc,3)), paste0("Microbiota, AUC=", round(roc2$auc,3)))) +
  theme(axis.text = element_text(size=11, face="bold")) + 
  theme(axis.title = element_text(size=11,face = "bold")) +
  labs(x="Specificity", y="Sensitivity")  +
  theme(
    legend.position = c(0.98, 0.1),
    legend.justification = c("right", "bottom"),
    legend.box.just = "right") + 
  theme(legend.title=element_blank(),
        legend.text = element_text(size=9, face="bold")) +
  theme(plot.title = element_text(face = "bold", size=14)) +
  ggtitle("ROC Curve")

print(q3)

#### 10.2) To visualize all decision trees 

#k is the index of the tree to be plotted. 
reprtree:::plot.getTree(classifier3, k=300)

######################
#### SUMMARY PLOT ####
######################

summary_plot <- ggarrange(ggarrange(p1,p2, ncol=2, labels=c("A","B")), 
                      ggarrange(p3, p4, p5, ncol = 3, labels = c("C", "D", "E"), widths = c(1,1,1)), 
                      ggarrange(q1,q2,q3,ncol=3, labels=c("F", "G", "H")),
                      nrow = 3, heights = c(5,3,4))         
ggsave(file="Summary_plot.png", summary_plot, height=12, width=15, units="in")
