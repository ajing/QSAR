library(corrplot)
library(caret)
#corrplot: the library to compute correlation matrix.

qsar <- read.table("QSAR data2.csv", header = TRUE, sep = ",")
qsar_no_na <- qsar[, colSums(!is.na(qsar)) == nrow(qsar)]
#read the tab file using the read table function.

qsar.scale<- scale(qsar_no_na,center=TRUE,scale=TRUE)
nzv <- nearZeroVar(qsar.scale)

qsar.scale_no_na <- qsar.scale[, colSums(!is.na(qsar.scale)) == nrow(qsar.scale)]
qsar.scale_no_na <- qsar.scale_no_na[,-c(1:4)]
# remove near zero variance variable also
nzv <- nearZeroVar(qsar.scale_no_na)
qsar.scale_no_na <- qsar.scale_no_na[,-nzv]

cor_qsar <- cor(qsar.scale_no_na, use = "complete.obs")
#compute the correlation matrix

cor_qsar <- cor(qsar_log, use = "complete.obs")
#compute the correlation matrix

highlyCorDescr <- findCorrelation(cor_qsar, cutoff = .70)
filteredDescr <- qsar.scale_no_na[,-highlyCorDescr]
descrCor2 <- cor(cbind(qsar.scale[,1:4], filteredDescr))
summary(descrCor2[upper.tri(descrCor2)])
#remove highly correlated features


pdf("cor.pdf")
corrplot(descrCor2, order = "hclust", tl.cex = 0.5)
dev.off()
#visualize the matrix, clustering features by correlation index.f


################## For those descriptor with high correlation  to Ct and Cp
qsar.scale_no_na <- qsar.scale[, colSums(is.na(qsar.scale)) != nrow(qsar.scale)]
cor_qsar <- cor(qsar.scale_no_na, use = "complete.obs")
col_desc <- rownames(cor_qsar)[rowSums(abs(cor_qsar[,c("Cp1", "Cp6", "Ct1", "Cp6", "Ratio1", "Ratio6")]) > 0.4) > 1]

pdf("cor_sel.pdf")
corrplot(cor_qsar[col_desc, col_desc], order = "hclust", tl.cex = 0.7)
dev.off()


cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1]
      uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}

res1 <- cor.mtest(cor_qsar[col_desc, col_desc], 0.95)
#res2 <- cor.mtest(mtcars, 0.99)
## specialized the insignificant value according to the significant level
pdf("cor_sig.pdf")
corrplot.mixed(cor_qsar[col_desc, col_desc], p.mat = res1[[1]], sig.level = 0.2, order = "hclust", tl.cex = 0.7)
dev.off()


############## features selection ##############
library(Boruta)
important <- Boruta(Ratio6 ~ ., data=data.frame(subset(qsar.scale_no_na, select = c(a_nF, vsurf_DW13, vsurf_ID8, CASA..1, ASA_P, vsurf_HB3, density, vsurf_HB8, vsurf_HB7, Ratio6))),  doTrace = 2, ntree = 30)
plot(important)


important <- Boruta(Ratio6 ~ ., data=data.frame(subset(qsar.scale_no_na, select = c(a_nF, vsurf_DW13, vsurf_ID8, CASA..1, ASA_P, vsurf_HB3, density, vsurf_HB8, vsurf_HB7, Ratio6))),  doTrace = 2, ntree = 30)
plot(important)

important <- Boruta(Ratio6 ~ ., data=data.frame(cbind(subset(qsar.scale_no_na, select = grep("vsurf", colnames(qsar.scale_no_na))), Ratio6= qsar.scale_no_na[,"Ratio6"])),  doTrace = 2, ntree = 30)
plot(important, las=2)

important <- Boruta(Ratio1 ~ ., data=data.frame(qsar.scale_no_na), doTrace = 2, ntree = 30)
plot(important, las=2)

important <- Boruta(Ct1 ~ ., data=data.frame(qsar.scale[, colSums(!is.na(qsar.scale)) == nrow(qsar.scale)]), doTrace = 2, ntree = 30)
plot(important, las=2)


lm_model_int <- lm(Ratio1 ~ density:a_nF + density:vsurf_DW13 + density:vsurf_ID8 + density:CASA..1 + density:ASA_P +  density:vsurf_HB3 + density:vsurf_HB8 + density:vsurf_HB7, data = data.frame(qsar))
summary(lm_model_int)



################# new data ####################
qsar <- read.table("QSARalldata.csv", header = TRUE, sep = ",")
interested <- names(important$finalDecision[important$finalDecision != "Rejected"])
ln_model <- lm(Ratio1 ~. , data = subset(qsar, select = colnames(qsar) %in% c(interested, "Ratio1")))
summary(ln_model)

ln_model <- lm(Ratio1 ~. , data = data.frame(subset(qsar.scale, select = colnames(qsar.scale) %in% c("a_nH","BCUT_PEOE_0", "vsurf_CW6", "vsurf_HL2", "Ratio1"))))
summary(ln_model)

rfFit1 <- train(x = subset(qsar.scale, select = colnames(qsar.scale) %in% c("a_nH","BCUT_PEOE_0", "vsurf_CW6", "vsurf_HL2")), y = qsar.scale[, "Ratio1"],
                 method = "rf",
                 preProcess = c("center", "scale"),
                 tuneLength = 10,
                 trControl = trainControl(method = "LOOCV", repeats = 5))


gbmFit1 <- train(x = subset(qsar, select = colnames(qsar) %in% c("a_nH","BCUT_PEOE_0", "vsurf_CW6", "vsurf_HL2")), y = qsar[, "Ratio1"],
                method = "gbm",
                preProcess = c("center", "scale"),
                tuneLength = 10,
                trControl = trainControl(method = "LOOCV", repeats = 5))

# xgboost
grid <- expand.grid(nrounds = seq(1, 201, by = 25),
                    max_depth = 1:6,
                    eta = (1:4)/10)
xgbFit1 <- train(Ratio1 ~ ., data = subset(qsar.scale[,-nzv], select = -Ratio6),
                 method = modelInfo,
                 tuneGrid = grid,
                 trControl = trainControl(method = 'LOOCV'))

> max(xgbFit1$results$Rsquared)
[1] 0.1460117

xgbFit2 <- train(log(Ratio1) ~ ., data = subset(qsar[qsar$Ratio1 > 0, ], select = c(vsurf_DW13, vsurf_ID8, CASA..1, ASA_P, vsurf_HB3, density, vsurf_HB8, vsurf_HB7, Ratio1)),
                 method = modelInfo,
                 tuneGrid = grid,
                 trControl = trainControl(method = 'LOOCV'))

tmp <- lm(Ct1 ~ ., data = subset(qsar, select = c(Ct1, Cp1, vsurf_DW13, vsurf_ID8, CASA..1, ASA_P, vsurf_HB3, density, vsurf_HB8, vsurf_HB7)))

#xgbFinal <- xgboost(param=xgbFit1$finalModel$tuneValue, data = subset(qsar.scale[,-nzv], select = -c(Ratio1, Ratio6)), label = qsar$Ratio1, nrounds = xgbFit1$finalModel$tuneValue$nrounds)


ctrl <- gafsControl(functions = caretGA)
obj <- gafs(x = subset(qsar.scale, select = colnames(qsar.scale) %in% c("a_nH","BCUT_PEOE_0", "vsurf_CW6", "vsurf_HL2", "Ratio1")),
            y = qsar.scale[, "Ratio1"],
            iters = 100,
            gafsControl = ctrl,
            ## Now pass options to `train`
            method = "")


ctrl <- safsControl(functions = caretSA)
obj <- safs(x = subset(qsar.scale, select = -c(Ratio1, Ratio6)),
            y = qsar.scale[, "Ratio1"],
            iters = 100,
            safsControl = ctrl,
            ## Now pass options to `train`
            method = "gbm")







####################################################
##############  Take Log   ##################
####################################################

# the number of samples
n_sample = 35

# take log for Cp Ct
qsar_log = qsar[1:n_sample, ]
qsar_log[, c("Ct1", "Cp1", "Ct6", "Cp6")] = log(qsar[1:n_sample, c("Ct1", "Cp1", "Ct6", "Cp6")] + 1)
nzv <- nearZeroVar(qsar_log)
qsar_log <- qsar_log[,-nzv]

# keep predictor features
predictors <- qsar_log[, c("Ct1", "Cp1", "Ct6", "Cp6", "Ratio1", "Ratio6")]

cor_qsar_all <- cor(qsar_log, use = "complete.obs")
#compute the correlation matrix

qsar_log <- subset(qsar_log, select=-c(Ct1, Cp1, Ct6, Cp6, Ratio1, Ratio6))
qsar_log <- qsar_log[, colSums(!is.na(qsar_log)) == nrow(qsar_log)]

cor_qsar <- cor(qsar_log, use = "complete.obs")
#compute the correlation matrix

#corrplot(cor_qsar, order = "hclust", tl.cex = 0.5)

highlyCorDescr <- findCorrelation(cor_qsar, cutoff = .65)
filteredDescr <- qsar_log[,-highlyCorDescr]
descrCor2 <- cor(cbind(predictors, filteredDescr))
summary(descrCor2[upper.tri(descrCor2)])
#remove highly correlated features


pdf("cor.pdf")
corrplot(descrCor2, order = "hclust", tl.cex = 0.5)
dev.off()
#visualize the matrix, clustering features by correlation index.f



#1. log(Ct), log(Cp) and other qsar features 的matrix correlation (only training data)
# done

#2. modeling log(Ct) ~ log(Cp) + other variables (1hour and 6 hour, respectively) (only training data)

#2.1 for 1 hour
# lasso for feature selection
y = predictors[,c("Ct1")]
x = as.matrix(cbind(predictors[,c("Cp1")], qsar_log))

library(glmnet)
cvfit <- glmnet(x = scale(x), y = scale(y))
plot(cvfit)
selected <- coef(cvfit, s = 0.07)
selected <- names(selected[,1][selected[,1] > 0])[-1]
selected <- selected[!(selected %in% c("ast_violation_ext", "lip_don", "vsurf_CW3", "vsurf_W7"))]

tmp_lm <- lm(Ct1 ~., data = data.frame(cbind(predictors[,c("Ct1", "Cp1")], qsar_log[, selected])))
summary(tmp_lm)

#2.1 for 6 hour
# lasso for feature selection
y = predictors[,c("Ct6")]
x = as.matrix(cbind(predictors[,c("Cp6")], qsar_log))

library(glmnet)
cvfit <- glmnet(x = scale(x), y = scale(y))
plot(cvfit)
selected <- coef(cvfit, s = 0.1)
selected <- names(selected[,1][selected[,1] > 0])[-1]
selected <- selected[!(selected %in% c("ast_violation_ext", "lip_don", "vsurf_CW3", "vsurf_W7"))]

tmp_lm <- lm(Ct6 ~., data = data.frame(cbind(predictors[,c("Ct6", "Cp6")], qsar_log[, selected])))
summary(tmp_lm)



#3. predict Ct at 1, 6 hour based on Cp 1 and / or Cp6 + other variables using the model in the step 2. (only test data)

#3.1 for 1 hour
# lasso for feature selection
y = predictors[,c("Ct1")]
x = as.matrix(cbind(predictors[,c("Cp1", "Cp6")], qsar_log))

library(glmnet)
cvfit <- glmnet(x = scale(x), y = scale(y))
plot(cvfit)
selected <- coef(cvfit, s = 0.08)
selected <- names(selected[,1][selected[,1] > 0])[-1]
selected <- selected[!(selected %in% c("ast_violation_ext", "lip_don", "vsurf_CW3", "vsurf_W7"))]

tmp_lm <- lm(Ct1 ~., data = data.frame(cbind(predictors[,c("Ct1", "Cp1", "Cp6")], qsar_log[, selected])))
summary(tmp_lm)

#3.1 for 1 hour
# lasso for feature selection
y = predictors[,c("Ct6")]
x = as.matrix(cbind(predictors[,c("Cp1", "Cp6")], qsar_log))

library(glmnet)
cvfit <- glmnet(x = scale(x), y = scale(y))
plot(cvfit)
selected <- coef(cvfit, s = 0.08)
selected <- names(selected[,1][selected[,1] > 0])[-1]
selected <- selected[!(selected %in% c("ast_violation_ext", "lip_don", "vsurf_CW3", "vsurf_W7"))]

tmp_lm <- lm(Ct6 ~., data = data.frame(cbind(predictors[,c("Ct6", "Cp1", "Cp6")], qsar_log[, selected])))
summary(tmp_lm)




#4. The residual between predicted Ct and experimental Ct.
