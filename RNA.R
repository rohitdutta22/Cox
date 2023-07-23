rna <- read.csv("METABRIC_RNA_Mutation.csv")
View(rna)
dim(rna)

rna <- na.omit(rna)  # removing the NA values from the data frame
dim(rna)

attach(rna)

##### type of breast surgery #####

pos1 <- which(type_of_breast_surgery == "")

pos2 <- which(cancer_type_detailed == "")

pos3 <- which(cellularity == "")

pos4 <- which(er_status_measured_by_ihc == "")

pos5 <- which(er_status_measured_by_ihc == "")

pos6 <- which(tumor_other_histologic_subtype == "")

pos7 <- which(primary_tumor_laterality == "")

pos8 <- which(oncotree_code == "")

pos9 <- which(X3.gene_classifier_subtype == "")

pos10 <- unique(c(which(death_from_cancer == ""),which(death_from_cancer == "Died of Other Causes")))





remove_rows <- c(pos1, pos2, pos3, pos4, pos5, pos6, pos7, pos8, pos9, pos10)
remove_rows <- unique(remove_rows)

library(dplyr)
rna <- rna %>% filter(!row_number() %in% remove_rows)     
dim(rna)
attach(rna)


factor_to_numeric_columns_index <- c(3,5,6,8,10,11,13,14,15,17,18,19,23,26,28,31)

for(i in 1:length(factor_to_numeric_columns_index)){
rna[,factor_to_numeric_columns_index[i]] <- as.numeric(as.factor(rna[,factor_to_numeric_columns_index[i]]))
}
View(rna)

rna$overall_survival <- abs(overall_survival - 1)
View(rna)
attach(rna)


drop_col <- c("patient_id", "cancer_type", "death_from_cancer")
rna = rna[,!(names(rna) %in% drop_col)]
View(rna)
attach(rna)


#### load rna_cleaned ####


dim(rna)
library(survival)
attach(rna)

drop_col_X <- c("overall_survival", "overall_survival_months")
rna <- rna[,!(names(rna) %in% drop_col_X)]
rna <- cbind(overall_survival, overall_survival_months, rna)

full_rna <- rna



attach(full_rna)

drop_col_X <- c("overall_survival", "overall_survival_months")
X <- full_rna[,!(names(full_rna) %in% drop_col_X)]
#View(X)








is.data.frame(X)
dim(X)


for(i in 516:688){
  X[,i] <- as.numeric(as.factor(X[,i]))
}
#View(X)


for(i in 518:690){
  full_rna[,i] <- as.numeric(as.factor(full_rna[,i]))
}



######################### detecting multicollinearity ######################3
ind <- array(0)
for(i in 1:(dim(X))[2]){
  if(sum(X[,i]) == 854){
    print(paste(i))
    ind <- c(ind,i)
  }
}
ind <- ind[-1]



X_cor <- X[,-ind]





cor_matrix <- cor(X_cor)
library(corrplot)
corrplot(cor_matrix)

upper <- array(0, dim = c(682,682))
for(i in 1:dim(cor_matrix)[1]){
  for(j in 1:dim(cor_matrix)[2]){
    if(i >= j){
      upper[i,j] <- 0 
    }else{
      upper[i,j] <- cor_matrix[i,j]
    }
  }
}
index_cor <- which(upper > 0.5 & cor_matrix < 1, arr.ind = TRUE)
dim(index_cor)

predictor_names <- colnames(X_cor)
predictor_names[as.vector(index_cor[,1])]

cor_upper <- data.frame(predictor1 = predictor_names[as.vector(index_cor[,1])],
                        predictor2 = predictor_names[as.vector(index_cor[,2])])
View(cor_upper)





index_cor <- which(upper > -1 & cor_matrix < -0.5, arr.ind = TRUE)
dim(index_cor)

predictor_names <- colnames(X_cor)
predictor_names[as.vector(index_cor[,1])]

cor_lower <- data.frame(predictor1 = predictor_names[as.vector(index_cor[,1])],
                        predictor2 = predictor_names[as.vector(index_cor[,2])])
View(cor_lower)





###################################################################





set.seed(1234)
split_index <- sample(1:dim(full_rna)[1], floor((dim(full_rna)[1])*0.4), replace = FALSE)
length(split_index)
rna <- full_rna[split_index,]
test.rna <- full_rna[-split_index,]


attach(rna)

cor_matrix <- cor(rna)
library(corrplot)
corrplot(cor_matrix)




##################### Lasso ###########################
################## lasso ##############################
library(glmnet)
df <- rna

#set.seed(4321)
#split_index <- sample(1:dim(df)[1], floor((dim(df)[1])*0.8), replace = FALSE)

test.glm <- test.rna
train.glm <- rna

X_predictors <- as.matrix(train.glm[,3:690])

fit <- glmnet(X_predictors, Surv(train.glm$overall_survival_months,
                                 train.glm$overall_survival),
              family = "cox", alpha = 1)
fit

cv.fit <- cv.glmnet(X_predictors, Surv(train.glm$overall_survival_months,
                                       train.glm$overall_survival),
                    family = "cox", alpha = 1)

cv.fit$lambda.min  #optimum lambda  


plot(cv.fit)

plot(cv.fit, xlim = c(-3,-1) , ylim = c(1,50))

cv.fit$lambda.1se  # minimum se 
fit_lasso <- glmnet(X_predictors, Surv(train.glm$overall_survival_months,
                                       train.glm$overall_survival),
                    family = "cox", alpha = 1, lambda = cv.fit$lambda.min)

coef <- fit_lasso$beta

nonzero_coef <- which(coef[,1] != 0)
features <- rownames(coef)[nonzero_coef]
features
index_lasso <- as.vector(which(coef[,1] != 0)) + 2

fit.lasso.final <- coxph(Surv(overall_survival_months,
                              overall_survival)~., train.glm[,c(1,2,index_lasso)])
summary(fit.lasso.final)

df_test <- test.glm[,c(1,2,index_lasso)]
(concordance(fit.lasso.final, newdata = df_test))$concordance
predicted_probs <- predict(fit.lasso.final, newdata = test.glm[,c(1,2,index_lasso)], type = "survival", se.fit = TRUE)

mean(predicted_probs$se.fit)




############################################################3
concordance_vec_lasso <- numeric(50)
nonzero_coef <- array(0)
index_error <- array(0)
error_count <- 0


set.seed(1234)
for(i in 1:50){
  print(paste("starting",i))
  
  split_index <- sample(1:dim(full_rna)[1], floor((dim(full_rna)[1])*0.4), replace = FALSE)
  rna <- full_rna[split_index,]
  test.rna <- full_rna[-split_index,]
  
  
  test.glm <- test.rna
  train.glm <- rna
  
  
  cv.fit <- cv.glmnet(X_predictors, Surv(train.glm$overall_survival_months,
                                         train.glm$overall_survival),
                      family = "cox", alpha = 1)
  
  fit_lasso <- glmnet(X_predictors, Surv(train.glm$overall_survival_months,
                                         train.glm$overall_survival),
                      family = "cox", alpha = 1, lambda = cv.fit$lambda.min)
  
  coef <- fit_lasso$beta
  
  nonzero_coef <- c(nonzero_coef,as.vector(which(coef[,1] != 0)))
  #features <- rownames(coef)[nonzero_coef]
  #features
  
  #features.coef <- coef[nonzero_coef]
  #names(features.coef) <- features
  #features.coef
  
  if(length(as.vector(which(coef[,1] != 0))) > 0){
    
    index_lasso <- as.vector(which(coef[,1] != 0)) + 2
    
    fit.lasso.final <- coxph(Surv(overall_survival_months,
                                  overall_survival)~., train.glm[,c(1,2,index_lasso)])
    
    concordance_vec_lasso[i] <- (concordance(fit.lasso.final, newdata = test.glm[,c(1,2,index_lasso)]))$concordance
  }else{
    
    fit.lasso.final <- coxph(Surv(overall_survival_months,
                                  overall_survival)~., train.glm)
    
    concordance_vec_lasso[i] <- (concordance(fit.lasso.final, newdata = test.glm))$concordance
    print("error may have occured")
    index_error <- c(index_error,i) 
    error_count <- error_count + 1
  }
  
}




mean(concordance_vec_lasso)
nonzero_coef <- nonzero_coef[-1]

tab <- table(nonzero_coef)
sort(tab, decreasing = TRUE)
index_error <- index_error[-1]
error_count

mean(concordance_vec_lasso[-index_error])



predictor_index <- as.numeric(rownames(sort(tab, decreasing = TRUE))) + 2
predictors <- (colnames(full_rna))[predictor_index]
predictor_count <- as.vector(sort(tab, decreasing = TRUE))
pred_data <- data.frame(predictor_index, predictors,
                        predictor_count)
save(pred_data, file = "pred_data_rna.Rdata")


#set.seed(1234)
#split_index <- sample(1:dim(df)[1], floor((dim(df)[1])*0.8), replace = FALSE)

test.glm <- test.rna
train.glm <- rna

concordance_vec_check <- numeric(length(predictor_index))
for(i in 1:length(predictor_index)){
  print(paste("starting",i))
  fit.lasso.check <- coxph(Surv(overall_survival_months,
                                overall_survival)~., train.glm[,c(1,2,predictor_index[1:i])])
  print(summary(fit.lasso.check))
  concordance_vec_check[i] <- (concordance(fit.lasso.check, newdata = test.glm[,c(1,2,predictor_index[1:i])]))$concordance
}

length(predictor_index)
length(concordance_vec_check)
plot(1:267,concordance_vec_check, type = "l",
     xlab = "Number of Predictors",
     ylab = "Value of Concordance")
points(which.max(concordance_vec_check),
       concordance_vec_check[which.max(concordance_vec_check)],
       col = "red", pch = 19)
concordance_vec_check[which.max(concordance_vec_check)]
which.max(concordance_vec_check)



