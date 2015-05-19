library(corrplot)
library(caret)
#corrplot: the library to compute correlation matrix.

qsar <- read.table("QSAR data2.csv", header = TRUE, sep = ",")
qsar_no_na <- qsar[, colSums(!is.na(qsar)) == nrow(qsar)]
#read the tab file using the read table function.

qsar.scale<- scale(qsar_no_na,center=TRUE,scale=TRUE);
qsar.scale_no_na <- qsar.scale[, colSums(is.na(qsar.scale)) != nrow(qsar.scale)]
qsar.scale_no_na <- qsar.scale_no_na[,-c(1:4)]
#scale all the features (from feature 2 bacause feature 1 is the predictor output)

cor_qsar <- cor(qsar.scale_no_na, use = "complete.obs")
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


lm_model_int <- lm(Ratio1 ~ density:a_nF + density:vsurf_DW13 + density:vsurf_ID8 + density:CASA..1 + density:ASA_P +  density:vsurf_HB3 + density:vsurf_HB8 + density:vsurf_HB7, data = data.frame(qsar))
summary(lm_model_int)
