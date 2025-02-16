#=======================================================


#=======================================================

rm(list=ls())
library(survival)         
library(Hmisc)            
library(dplyr)            
library(tidyverse)        
library(CoxBoost)         
library(randomForestSRC)  
library(gbm)               
library(caret)            
library(survivalsvm)     
library(glmnet)
library(superpc)
library(pheatmap)
library(timeROC)       
library(survivalROC)

#=======================================================

#=======================================================

train_data <- read.csv("traindata.csv", header = TRUE, row.names = 1)
test_data <- read.csv("testdata.csv", header = TRUE, row.names = 1)

#=======================================================

#=======================================================

model_names <- c("SRF_100", "SRF_500", "CoxPH", "CoxGB_100", "CoxGB_500", 
                 "GBM_100", "GBM_500", "sCox", "LASSO", "Ridge", 
                 "ENet_0.1", "ENet_0.5", "ENet_0.9")

# Create dataframe to store the scores
scores_df <- data.frame(matrix(NA, nrow = nrow(test_data), ncol = length(model_names)))
colnames(scores_df) <- model_names


#=======================================================
# RandomForestSRC with 100 trees and 500 trees
#=======================================================
set.seed(1234)

# rfSRC with 100 trees
model_RF_100 <- rfsrc(Surv(OSTime, OS) ~ ., data = train_data, ntree = 100)
scores_df[, "SRF_100"] <- predict(model_RF_100, test_data)$predicted

# rfSRC with 500 trees
model_RF_500 <- rfsrc(Surv(OSTime, OS) ~ ., data = train_data, ntree = 500)
scores_df[, "SRF_500"] <- predict(model_RF_500, test_data)$predicted

#=======================================================
# CoxPH (No hyperparameters to tune, just use default)
#=======================================================
fit_cox <- coxph(Surv(OSTime, OS) ~ ., data = train_data)
scores_df[, "CoxPH"] <- as.numeric(predict(fit_cox, newdata = test_data, type = "lp"))

#=======================================================
# CoxBoost with 100 and 200 boosting steps
#=======================================================
fit_CoxBoost_100 <- CoxBoost(time = train_data$OSTime, status = train_data$OS, x = as.matrix(train_data[, -c(1, 2)]), stepno = 100)
scores_df[, "CoxGB_100"] <- as.numeric(predict(fit_CoxBoost_100, newdata = as.matrix(test_data[, -c(1, 2)]), type = "lp"))

fit_CoxBoost_500 <- CoxBoost(time = train_data$OSTime, status = train_data$OS, x = as.matrix(train_data[, -c(1, 2)]), stepno = 500)
scores_df[, "CoxGB_500"] <- as.numeric(predict(fit_CoxBoost_500, newdata = as.matrix(test_data[, -c(1, 2)]), type = "lp"))

#=======================================================
# GBM with 100 and 500 trees
#=======================================================
fit_gbm_100 <- gbm(Surv(OSTime, OS) ~ ., data = train_data, distribution = 'coxph', n.trees = 100)
scores_df[, "GBM_100"] <- predict(fit_gbm_100, newdata = test_data, type = "link", n.trees = 100)

fit_gbm_500 <- gbm(Surv(OSTime, OS) ~ ., data = train_data, distribution = 'coxph', n.trees = 500)
scores_df[, "GBM_500"] <- predict(fit_gbm_500, newdata = test_data, type = "link", n.trees = 500)

#=======================================================
# Stepwise CoxPH (No hyperparameters to tune, just use default)
#=======================================================
full_cox_model <- coxph(Surv(OSTime, OS) ~ ., data = train_data)
step_cox_model <- step(full_cox_model, direction = "both", trace = FALSE)
scores_df[, "sCox"] <- as.numeric(predict(step_cox_model, newdata = test_data, type = "lp"))

#=======================================================
# Lasso, Ridge, and Elastic Net (Cross-validated lambda)
#=======================================================
x_train <- model.matrix(Surv(OSTime, OS) ~ ., data = train_data)[, -1]
y_train <- Surv(train_data$OSTime, train_data$OS)
x_test <- model.matrix(Surv(OSTime, OS) ~ ., data = test_data)[, -1]

# Lasso (alpha = 1)
fit_lasso <- cv.glmnet(x_train, y_train, family = "cox", alpha = 1)
scores_df[, "LASSO"] <- as.numeric(predict(fit_lasso, newx = x_test, s = "lambda.min", type = "link"))

# Ridge (alpha = 0)
fit_ridge <- cv.glmnet(x_train, y_train, family = "cox", alpha = 0)
scores_df[, "Ridge"] <- as.numeric(predict(fit_ridge, newx = x_test, s = "lambda.min", type = "link"))

# Elastic Net with alpha = 0.1
fit_enet_0.1 <- cv.glmnet(x_train, y_train, family = "cox", alpha = 0.1)
scores_df[, "ENet_0.1"] <- as.numeric(predict(fit_enet_0.1, newx = x_test, s = "lambda.min", type = "link"))

# Elastic Net with alpha = 0.5
fit_enet_0.5 <- cv.glmnet(x_train, y_train, family = "cox", alpha = 0.5)
scores_df[, "ENet_0.5"] <- as.numeric(predict(fit_enet_0.5, newx = x_test, s = "lambda.min", type = "link"))

# Elastic Net with alpha = 0.9
fit_enet_0.9 <- cv.glmnet(x_train, y_train, family = "cox", alpha = 0.9)
scores_df[, "ENet_0.9"] <- as.numeric(predict(fit_enet_0.9, newx = x_test, s = "lambda.min", type = "link"))

#=======================================================

#=======================================================

scores_df <- as.data.frame(scale(scores_df))
rownames(scores_df) <- rownames(test_data)

column_names <- names(scores_df)

for (i in 1:(length(column_names) - 1)) {
  for (j in (i + 1):length(column_names)) {
    new_column_name <- paste(column_names[i], column_names[j], sep="_plus_")
    scores_df[[new_column_name]] <- scores_df[[column_names[i]]] + scores_df[[column_names[j]]]
  }
}


scores_df$OSTime <- test_data$OSTime
scores_df$OS <- test_data$OS

#=======================================================

#=======================================================

column_names <- names(scores_df)
filtered_column_names <- column_names[!column_names %in% c("OS", "OSTime")]


model_timeAUCs <- numeric(length(filtered_column_names))

for (i in 1:length(filtered_column_names)) {
  model_name <- filtered_column_names[i]
  risk_scores <- scores_df[, model_name]
  
  
  roc_result <- survivalROC(Stime = scores_df$OSTime,
                            status = scores_df$OS,
                            marker = risk_scores,
                            predict.time = 5,
                            method = "KM")
  
  timeauc <- roc_result$AUC
  print(timeauc)
  model_timeAUCs[i] <- timeauc
  print(i)
}

timeAUC_df <- data.frame(Model = filtered_column_names, timeAUC = model_timeAUCs)

timeAUC_df

write.csv(timeAUC_df, file = "timeAUC_test.csv")


#=======================================================





#=======================================================






#=======================================================



rm(list=ls())


selected_data <- read.csv("traindata.csv", header = T, row.names = 1)



set.seed(1234)
random_order <- sample(nrow(selected_data)) 
shuffled_data <- selected_data[random_order, ]
print(shuffled_data[1:4, 1:4])
selected_data <- shuffled_data



K <- 5

model_names <-  c("SRF_100", "SRF_500", "CoxPH", "CoxGB_100", "CoxGB_500", 
                  "GBM_100", "GBM_500", "sCox", "LASSO", "Ridge", 
                  "ENet_0.1", "ENet_0.5", "ENet_0.9")


#=======================================================

#=======================================================


auc_list <- list()

for (fold in 1:K) {
  
  print(fold)
  fold_size <- nrow(selected_data) / K
  test_indices <- ((fold - 1) * fold_size + 1):(fold * fold_size)
  train_indices <- setdiff(1:nrow(selected_data), test_indices)
  train_data <- selected_data[train_indices, ]
  test_data <- selected_data[test_indices, ]
  

  model_RF_100 <- rfsrc(Surv(OSTime, OS) ~ ., data = train_data, ntree = 100)
  risk_scores_RF_100 <- predict(model_RF_100, test_data)$predicted
  model_RF_500 <- rfsrc(Surv(OSTime, OS) ~ ., data = train_data, ntree = 500)
  risk_scores_RF_500 <- predict(model_RF_500, test_data)$predicted
  print("rf model construction completed")
  current_time <- Sys.time()
  formatted_time <- format(current_time, "%Y-%m-%d %H:%M:%S")
  print(formatted_time)
  
  

  fit_cox <- coxph(Surv(OSTime, OS) ~ ., data = train_data)
  risk_scores_cox <- predict(fit_cox, newdata = test_data[, -c(1, 2)], type = "lp")
  risk_scores_cox <- as.numeric(risk_scores_cox)
  print("cox model construction completed")
  current_time <- Sys.time()
  formatted_time <- format(current_time, "%Y-%m-%d %H:%M:%S")
  print(formatted_time)
  
  
  
  fit_CoxBoost_100 <- CoxBoost(train_data[, 'OSTime'], train_data[, 'OS'], as.matrix(train_data[, -c(1, 2)]),
                               stepno = 100)
  risk_scores_CoxBoost_100 <- predict(fit_CoxBoost_100, newdata = test_data[, -c(1, 2)], type = "lp")
  risk_scores_CoxBoost_100 <- as.numeric(risk_scores_CoxBoost_100)
  fit_CoxBoost_500 <- CoxBoost(train_data[, 'OSTime'], train_data[, 'OS'], as.matrix(train_data[, -c(1, 2)]),
                               stepno = 500)
  risk_scores_CoxBoost_500 <- predict(fit_CoxBoost_500, newdata = test_data[, -c(1, 2)], type = "lp")
  risk_scores_CoxBoost_500 <- as.numeric(risk_scores_CoxBoost_500)
  print("CoxBoost model construction completed")
  current_time <- Sys.time()
  formatted_time <- format(current_time, "%Y-%m-%d %H:%M:%S")
  print(formatted_time)
  
  
  fit_gbm_100 <- gbm(formula = Surv(OSTime, OS) ~ ., data = train_data, distribution = 'coxph')
  risk_scores_gbm_100 <- predict(fit_gbm_100, newdata = test_data)
  fit_gbm_500 <- gbm(formula = Surv(OSTime, OS) ~ ., data = train_data, distribution = 'coxph')
  risk_scores_gbm_500 <- predict(fit_gbm_500, newdata = test_data)
  print("gbm model construction completed")
  current_time <- Sys.time()
  formatted_time <- format(current_time, "%Y-%m-%d %H:%M:%S")
  print(formatted_time)
  
  
  
  full_cox_model <- coxph(Surv(OSTime, OS) ~ ., data = train_data)
  step_cox_model <- step(full_cox_model, direction = "both", trace = FALSE)
  risk_scores_stepcox <- predict(step_cox_model, newdata = test_data[, -c(1, 2)], type = "lp")
  risk_scores_stepcox <- as.numeric(risk_scores_stepcox)
  print("stepwise Cox model construction completed")
  current_time <- Sys.time()
  formatted_time <- format(current_time, "%Y-%m-%d %H:%M:%S")
  print(formatted_time)
  

  x_train <- model.matrix(Surv(OSTime, OS) ~ ., data = train_data)[, -1]
  y_train <- Surv(train_data$OSTime, train_data$OS)
  x_test <- model.matrix(Surv(OSTime, OS) ~ ., data = test_data)[, -1]
  
  fit_lasso <- cv.glmnet(x_train, y_train, family = "cox", alpha = 1)
  risk_scores_lasso <- predict(fit_lasso, newx = x_test, s = "lambda.min", type = "link")
  risk_scores_lasso <- as.numeric(risk_scores_lasso)
  print("lasso model construction completed")
  current_time <- Sys.time()
  formatted_time <- format(current_time, "%Y-%m-%d %H:%M:%S")
  print(formatted_time)
  

  fit_ridge <- cv.glmnet(x_train, y_train, family = "cox", alpha = 0)
  risk_scores_ridge <- predict(fit_ridge, newx = x_test, s = "lambda.min", type = "link")
  risk_scores_ridge <- as.numeric(risk_scores_ridge)
  print("ridge model construction completed")
  current_time <- Sys.time()
  formatted_time <- format(current_time, "%Y-%m-%d %H:%M:%S")
  print(formatted_time)
  
  
  
  fit_enet_0.1 <- cv.glmnet(x_train, y_train, family = "cox", alpha = 0.1)
  risk_scores_enet_0.1 <- predict(fit_enet_0.1, newx = x_test, s = "lambda.min", type = "link")
  risk_scores_enet_0.1 <- as.numeric(risk_scores_enet_0.1)
  
  fit_enet_0.5 <- cv.glmnet(x_train, y_train, family = "cox", alpha = 0.5)
  risk_scores_enet_0.5 <- predict(fit_enet_0.5, newx = x_test, s = "lambda.min", type = "link")
  risk_scores_enet_0.5 <- as.numeric(risk_scores_enet_0.5)
  
  fit_enet_0.9 <- cv.glmnet(x_train, y_train, family = "cox", alpha = 0.9)
  risk_scores_enet_0.9 <- predict(fit_enet_0.9, newx = x_test, s = "lambda.min", type = "link")
  risk_scores_enet_0.9 <- as.numeric(risk_scores_enet_0.9)
  
  
  print("Elastic Net model construction completed")
  current_time <- Sys.time()
  formatted_time <- format(current_time, "%Y-%m-%d %H:%M:%S")
  print(formatted_time)


  
  scores_df <- data.frame(matrix(NA, nrow = nrow(test_data), ncol = length(model_names)))
  colnames(scores_df) <- model_names
  
  scores_df[, "SRF_100"] <- risk_scores_RF_100
  scores_df[, "SRF_500"] <- risk_scores_RF_500
  scores_df[, "CoxPH"] <- risk_scores_cox
  scores_df[, "CoxGB_100"] <- risk_scores_CoxBoost_100
  scores_df[, "CoxGB_500"] <- risk_scores_CoxBoost_500
  
  scores_df[, "GBM_100"] <- risk_scores_gbm_100
  scores_df[, "GBM_500"] <- risk_scores_gbm_500
  scores_df[, "sCox"] <- risk_scores_stepcox
  scores_df[, "LASSO"] <- risk_scores_lasso
  scores_df[, "Ridge"] <- risk_scores_ridge
  
  scores_df[, "ENet_0.1"] <- risk_scores_enet_0.1
  scores_df[, "ENet_0.5"] <- risk_scores_enet_0.5
  scores_df[, "ENet_0.9"] <- risk_scores_enet_0.9
  
  #=======================================================
  
  #=======================================================
  
  scores_df <- as.data.frame(scale(scores_df))
  rownames(scores_df) <- rownames(test_data)
  
  column_names <- names(scores_df)
  
  for (i in 1:(length(column_names) - 1)) {
    for (j in (i + 1):length(column_names)) {
      new_column_name <- paste(column_names[i], column_names[j], sep="|")
      scores_df[[new_column_name]] <- scores_df[[column_names[i]]] + scores_df[[column_names[j]]]
    }
  }
  
  
  scores_df$OSTime <- test_data$OSTime
  scores_df$OS <- test_data$OS
  
  
  
  column_names <- names(scores_df)
  filtered_column_names <- column_names[!column_names %in% c("OS", "OSTime")]
  
  #=======================================================
  
  #=======================================================
  
  
  model_timeAUCs <- numeric(length(filtered_column_names))
  
  for (i in 1:length(filtered_column_names)) {
    model_name <- filtered_column_names[i]
    risk_scores <- scores_df[, model_name]
    
    
    
    roc_result <- survivalROC(Stime = scores_df$OSTime,
                              status = scores_df$OS,
                              marker = risk_scores,
                              predict.time = 5,
                              method = "KM")
    timeauc <- roc_result$AUC
    model_timeAUCs[i] <- timeauc
    print(i)
    print(timeauc)
  }
  
  timeAUC_df <- data.frame(Model = filtered_column_names, timeAUC = model_timeAUCs)
  
  auc_list[[fold]] <- timeAUC_df
  
}

auc_list


#=======================================================

#=======================================================


auc_columns <- lapply(auc_list, function(df) df$timeAUC)
model_names <- auc_list[[1]]$Model
combined_auc_df <- data.frame(Model = model_names)
combined_auc_df <- cbind(combined_auc_df, do.call(cbind, auc_columns))
print(combined_auc_df)
combined_auc_df$Model <- gsub("_", "", combined_auc_df$Model)


colnames(combined_auc_df) <- paste0('ValidationSet', colnames(combined_auc_df))


testResult <- read.csv("timeAUC_test.csv", header = T, row.names = 1)
rownames(testResult) == rownames(combined_auc_df)
combined_auc_df$testingSet <- testResult$timeAUC
rownames(combined_auc_df) <- combined_auc_df$ValidationSetModel
combined_auc_df$ValidationSetModel <- NULL
combined_auc_df$meanValue <- rowMeans(combined_auc_df, na.rm = TRUE)
combined_auc_df <- combined_auc_df %>%
  dplyr::arrange(desc(meanValue))

head(combined_auc_df)

# combined_auc_df <- round(combined_auc_df, 3)
# rownames(combined_auc_df) <- paste0(rownames(combined_auc_df), "(AUC:", combined_auc_df$meanValue, ")")


#=======================================================

#=======================================================


minVt <- min(combined_auc_df)
maxVt <- max(combined_auc_df)

colors <- colorRampPalette(c("#F6F5F2", "#B1C381"))(50)

breaks <- seq(minVt, maxVt , length.out=50) 


p <- pheatmap(as.matrix(combined_auc_df), 
              scale="none",
              display_numbers=TRUE,
              color=colors,
              breaks=breaks,
              cluster_rows=FALSE,  
              cluster_cols=FALSE,  
              width=8,  
              height=6, 
              main="AUC in year5",
              legend = FALSE,
              fontsize=7,  
              border_color=NA, 
              number_color = "black",  
              angle_col=45,
              number_format = "%.3f",
              na_col="gray")  

print(p)

pdf(file = "auc5.pdf", height = 12, width = 8)
print(p)
dev.off()


