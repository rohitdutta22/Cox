library(survival)
library(penalized)


#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("RLassoCox")



data("nki70")
df <- nki70

df$Grade <- as.numeric(df$Grade)

df$Diam <- as.numeric(as.factor(df$Diam))

df$N <- as.numeric(as.factor(df$N))

df$ER <- as.numeric(as.factor(df$ER))

attach(df)

head(df)


i = 1
set.seed(1234)

concordance.vec <- numeric(100)
se.vec <- numeric(100)
for(i in 1:100){
  print(paste("starting :",i))
  split_index <- sample(1:dim(df)[1], floor((dim(df)[1])*0.8), replace = FALSE)
  
  test <- df[-split_index,]
  train <- df[split_index,]
  
  #head(test)
  #head(train)
  
  fit <- coxph(Surv(time, event)~., train) 
  
  model.results <- summary(fit)
  
  y <- model.results$coefficients
  
  index <- as.vector(which(y[,5] <= 0.05)) + 2
  
  train.new <- train[,c(1,2,index)]
  test.new <- test[,c(1,2,index)]
  
  fit.new <- coxph(Surv(time, event)~., train.new)
  #print(summary(fit.new))
  
 # con <- (concordance(fit.new, newdata = test.new))$concordance
  
  concordance.vec[i] <- (concordance(fit.new, newdata = test.new))$concordance
  
  predicted_probs <- predict(fit.new, newdata = test.new, type = "survival", se.fit = TRUE)
  
  se.vec[i] <- mean(predicted_probs$se.fit)
  
  i = i + 1
}
mean(concordance.vec[1:47])
se.vec[1:47]
  
  



  
################## Stepwise selection in coxph ##########################
library(My.stepwise)
set.seed(4321)
split_index <- sample(1:dim(df)[1], floor((dim(df)[1])*0.8), replace = FALSE)
  
test <- df[-split_index,]
train <- df[split_index,]

my.variable.list <- (colnames(train))[3:77]
mod <- My.stepwise.coxph(Time = "time", Status = "event", variable.list = my.variable.list,
                        data = train)

fit.stepwise <- coxph(Surv(time, event)~ KNTC2 + ESM1 + ER +
                   Contig32125_RC + RUNDC1 + COL4A2 + PITRM1 + 
                   ORC6L + GPR126 + SLC2A3 + STK32B +
                   MMP9 + GPR180 + RTN4RL1 + RAB6B + SERF1A +
                   QSCN6L1 + Age + SCUBE2 + MELK + IGFBP5.1, train)

summary(fit.stepwise)

(concordance(fit.stepwise, newdata = test))$concordance

predicted_probs <- predict(fit.stepwise, newdata = test, type = "survival", se.fit = TRUE)

mean(predicted_probs$se.fit)






################## lasso ##############################
library(glmnet)

set.seed(1234)
split_index <- sample(1:dim(df)[1], floor((dim(df)[1])*0.8), replace = FALSE)

test.glm <- df[-split_index,]
train.glm <- df[split_index,]

X_predictors <- as.matrix(train.glm[,3:77])





##################################### looking for multicolinearity######################3
cor_matrix <- cor(X_predictors)
library(corrplot)
corrplot(cor_matrix)
dim(cor_matrix)

upper <- array(0, dim = c(75,75))
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

predictor_names <- colnames(X_predictors)
predictor_names[as.vector(index_cor[,1])]

cor_upper <- data.frame(predictor1 = predictor_names[as.vector(index_cor[,1])],
                        predictor2 = predictor_names[as.vector(index_cor[,2])])
View(cor_upper)





index_cor <- which(upper > -1 & cor_matrix < -0.5, arr.ind = TRUE)
dim(index_cor)

predictor_names <- colnames(X_predictors)
predictor_names[as.vector(index_cor[,1])]

cor_lower <- data.frame(predictor1 = predictor_names[as.vector(index_cor[,1])],
                        predictor2 = predictor_names[as.vector(index_cor[,2])])
View(cor_lower)


################################################################################








fit <- glmnet(X_predictors, Surv(train.glm$time, train.glm$event), family = "cox", alpha = 1)
fit

cv.fit <- cv.glmnet(X_predictors, Surv(train.glm$time, train.glm$event),
                    family = "cox", alpha = 1)

cv.fit$lambda.min  #optimum lambda  0.1438282, 0.05552002


plot(cv.fit)

plot(cv.fit, xlim = c(-3,-1) , ylim = c(1,50))

cv.fit$lambda.1se  # minimum se  0.2086702, 0.2042237
fit_lasso <- glmnet(X_predictors, Surv(train.glm$time, train.glm$event),
                    family = "cox", alpha = 1, lambda = cv.fit$lambda.min)

coef <- fit_lasso$beta

nonzero_coef <- which(coef[,1] != 0)
features <- rownames(coef)[nonzero_coef]
features
index_lasso <- as.vector(which(coef[,1] != 0)) + 2

fit.lasso.final <- coxph(Surv(time, event)~., train.glm[,c(1,2,index_lasso)])
summary(fit.lasso.final)

(concordance(fit.lasso.final, newdata = test.glm[,c(1,2,index_lasso)]))$concordance
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
split_index <- sample(1:dim(df)[1], floor((dim(df)[1])*0.8), replace = FALSE)

test.glm <- df[-split_index,]
train.glm <- df[split_index,]

cv.fit <- cv.glmnet(X_predictors, Surv(train.glm$time, train.glm$event),
                    family = "cox", alpha = 1)

fit_lasso <- glmnet(X_predictors, Surv(train.glm$time, train.glm$event),
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

  fit.lasso.final <- coxph(Surv(time, event)~., train.glm[,c(1,2,index_lasso)])
  
  concordance_vec_lasso[i] <- (concordance(fit.lasso.final, newdata = test.glm[,c(1,2,index_lasso)]))$concordance
}else{
  
  fit.lasso.final <- coxph(Surv(time, event)~., train.glm)
  
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
predictors <- (colnames(df))[predictor_index]
predictor_count <- as.vector(sort(tab, decreasing = TRUE))
pred_data <- data.frame(predictor_index, predictors,
                        predictor_count)
save(pred_data, file = "pred_data.Rdata")

set.seed(1234)
split_index <- sample(1:dim(df)[1], floor((dim(df)[1])*0.8), replace = FALSE)

test.glm <- df[-split_index,]
train.glm <- df[split_index,]

concordance_vec_check <- numeric(length(predictor_index))
for(i in 1:length(predictor_index)){
  print(paste("starting",i))
  fit.lasso.check <- coxph(Surv(time, event)~., train.glm[,c(1,2,predictor_index[1:i])])
  print(summary(fit.lasso.check))
  concordance_vec_check[i] <- (concordance(fit.lasso.check, newdata = test.glm[,c(1,2,predictor_index[1:i])]))$concordance
}

plot(1:55,concordance_vec_check, type = "l",
     xlab = "Number of Predictors",
     ylab = "Value of Concordance")
points(which.max(concordance_vec_check),
       concordance_vec_check[which.max(concordance_vec_check)],
       col = "red", pch = 19)
concordance_vec_check[which.max(concordance_vec_check)]



