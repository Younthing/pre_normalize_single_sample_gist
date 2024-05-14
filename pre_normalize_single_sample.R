##### 读取文件 #####
library(limma)
library(dplyr)
rm(list = ls())
exprs_126094 <- read.csv("./1-GEO下载/exprs_126094.csv", row.names = 1)
group_126094 <- read.csv("./1-GEO下载/group_126094.csv")
2
##### 是否log2（value+1）####
log2_transform <- function(matrix) { # 来自GEO官网
  print(paste0(c("最小值：", "最大值："), range(matrix)))
  qx <- as.numeric(quantile(matrix,
    c(0., 0.25, 0.5, 0.75, 0.99, 1.0),
    na.rm = TRUE
  ))
  log_c <- (qx[5] > 100) ||
    (qx[6] - qx[1] > 50 && qx[2] > 0)
  if (log_c) {
    matrix[which(matrix <= 0)] <- NaN
    matrix_log <- log2(matrix + 1)
    print("转换完成")
  } else {
    matrix_log <- matrix
    print("不需要转换")
  }
  return(matrix_log)
}
## log化
matrix_126094_log2 <- log2_transform(exprs_126094)
matrix_32323_log2 <- log2_transform(exprs_32323)

## 去批次----
library(sva)

## sva中batch指定批次，mod指定保护的差异
batch <- rep(1:10, 2)
group_126094$group <- factor(group_126094$group, levels = c("Normal", "CRC"))
mod <- model.matrix(~group, data = group_126094)
sva_batch_126094 <- sva::ComBat(exprs_126094, batch = batch, mod = mod)
limma_batch_126094 <- limma::removeBatchEffect(exprs_126094, batch = batch)

## normalization ----
## limma包normalizeBetweenArrays
## 先去批次再归一化
batch_126094_normal <- normalizeBetweenArrays(sva_batch_126094) %>%
  data.frame()
## 直接归一化-究极不推荐
matrix_126094_normal <- normalizeBetweenArrays(matrix_126094_log2) %>%
  data.frame()

## 可视化 ----
## 箱图
boxplot(matrix_126094_log2, outline = FALSE, las = 2) ## 看看样本间是否有批次
boxplot(sva_batch_126094, outline = FALSE, las = 2)
boxplot(batch_126094_normal, outline = FALSE, las = 2)

## 密度图
limma::plotDensities(matrix_126094_log2, legend = FALSE)
limma::plotDensities(sva_batch_126094, legend = FALSE)
limma::plotDensities(batch_126094_normal, legend = FALSE)

# PCA图  fviz_pca_ind
pca_before <- prcomp(t(matrix_126094_log2),
  scale = TRUE
)
pca_after <- prcomp(t(sva_batch_126094),
  scale = TRUE
)
pca_after_normalize <- prcomp(t(batch_126094_normal),
  scale = TRUE
)

group_col <- c("#5E86C1", "#f59292")
library(factoextra)
## 这是ggplot对象的
fviz_pca_ind(pca_before,
  title = "Before",
  label = "none",
  axes = c(1, 2),
  habillage = group_126094$group, # 按分组标注
  palette = group_col,
  addEllipses = TRUE, # 添加椭圆圈圈
  ellipse.level = 0.95 # 95％置信区间
)
## 仅去批次
fviz_pca_ind(pca_after,
  title = "After",
  label = "none",
  axes = c(1, 2),
  habillage = group_126094$group, # 按分组标注
  palette = group_col,
  addEllipses = TRUE, # 添加椭圆圈圈
  ellipse.level = 0.95 # 95％置信区间
)
## 去批次后再归一化的
fviz_pca_ind(pca_after_normalize,
  title = "After",
  label = "none",
  axes = c(1, 2),
  habillage = group_126094$group, # 按分组标注
  palette = group_col,
  addEllipses = TRUE, # 添加椭圆圈圈
  ellipse.level = 0.95 # 95％置信区间
)

cluster_data <- t(sva_batch_126094)
rownames(cluster_data) <- with(
  group_126094,
  paste(group, sample, sep = "-") # with
)
hc <- hclust(dist(cluster_data))
plot(hc)

cluster_data <- t(normalizeBetweenArrays(sva_batch_126094))
rownames(cluster_data) <- with(
  group_126094,
  paste(group, sample, sep = "-") # with
)
hc <- hclust(dist(cluster_data))
plot(hc)

## 样本相关性矩阵 ----
corr_matrix <- cor(normalizeBetweenArrays(limma_batch))
ggcorrplot::ggcorrplot(corr_matrix, lab = TRUE, lab_size = 2)

range(corr_matrix)