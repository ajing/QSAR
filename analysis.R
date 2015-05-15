library(corrplot)
library(caret)
#corrplot: the library to compute correlation matrix.

qsar <- read.table("QSAR data2.csv", header = TRUE, sep = ",")
qsar_no_na <- qsar[, colSums(is.na(qsar)) != nrow(qsar)]
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
col_desc <- rownames(cor_qsar)[rowSums(abs(cor_qsar[,c("Cp1", "Cp6", "Ct1", "Cp6")]) > 0.4) > 1]

pdf("cor_sel.pdf")
corrplot.mixed(cor_qsar[col_desc, col_desc], order = "hclust", tl.cex = 0.7)
dev.off()
