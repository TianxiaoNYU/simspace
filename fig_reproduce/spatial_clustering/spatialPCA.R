library(SpatialPCA)
library(ggplot2)
library(tidyr)
library(dplyr)

{
  rawcount <- read.csv(
    '../../fig_data/clustering/simspace_count.csv',
    header = T,
    row.names = 1)
  location <- read.csv(
    '../../fig_data/clustering/simspace_metadata.csv',
    header = T,
    row.names = 1)
  rawcount <- t(rawcount)
  location <- as.matrix(location)
  
  rownames(location) <- colnames(rawcount)
  new_location <- location[,c("row", "col")]
  
  real_rawcount <- read.csv('../../data/reference_count.csv',
                            header = T,
                            row.names = 1)
  real_meta <- read.csv('../../data/reference_metadata.csv',
                    header = T)
  real_location <- (real_meta[,c("x_centroid", "y_centroid")])
  # normalize the coordinates in real_location so that they sit between 0 and 100 for each column seperately
  real_location$x_centroid <- (real_location[,1] - min(real_location[,1])) / (max(real_location[,1]) - min(real_location[,1])) * 100
  real_location$y_centroid <- (real_location[,2] - min(real_location[,2])) / (max(real_location[,2]) - min(real_location[,2])) * 100
  colnames(real_location) <- c("row", "col")
  real_location$row <- real_location$row + 200
  rownames(real_location) <- paste0("real_", seq(1, nrow(real_location)))
  colnames(real_rawcount) <- rownames(real_location)
  real_rawcount <- as.matrix(real_rawcount)
  
  ## Load scCube data
  {
    scCube_rawcount <- read.csv('../../fig_data/clustering/scCube_count.csv', 
                         header = T,
                         row.names = 1)
    scCube_location <- read.csv('../../fig_data/clustering/scCube_metadata.csv',
                         header = T,
                         row.names = 1)
    scCube_location$point_x <- scCube_location$point_x + 100
    scCube_location <- as.matrix(scCube_location)
    
    rownames(scCube_location) <- colnames(scCube_rawcount)
    
    new_scCube_location <- scCube_location[,c("point_x", "point_y")]
    colnames(new_scCube_location) <- c("row", "col")
  }
  
  combined_rawcount <- cbind(rawcount, real_rawcount, scCube_rawcount)
  combined_location <- rbind(new_location, real_location, new_scCube_location)
  
  combined_rawcount <- as.matrix(combined_rawcount)
  combined_location <- as.matrix(combined_location)
  class(combined_location) <- "numeric"
}

## This might take a long time (up to several hours)
ST = CreateSpatialPCAObject(
  counts=combined_rawcount, 
  location=combined_location, 
  project = "SpatialPCA",
  gene.type="spatial",
  sparkversion="spark", 
  gene.number=100,
  numCores_spark = 8,
  customGenelist=NULL,
  min.loctions=1, 
  min.features=20)
ST = SpatialPCA_buildKernel(ST, kerneltype="gaussian", bandwidthtype="SJ")
ST = SpatialPCA_EstimateLoading(ST,fast=FALSE,SpatialPCnum=20)
ST = SpatialPCA_SpatialPCs(ST, fast=FALSE)

# saveRDS(ST, '/Users/zhaotianxiao/Library/CloudStorage/Dropbox/FenyoLab/Project/Spatialsim/reproduce/clustering/SpatialPCA/SpatialPCAres_0425.rds')
# ST <- readRDS('/Users/zhaotianxiao/Library/CloudStorage/Dropbox/FenyoLab/Project/Spatialsim/reproduce/clustering/SpatialPCA/SpatialPCAres_0425.rds')

combined_location <- as.data.frame(combined_location)
combined_location$Cluster <- c(location[,5], real_meta$Cluster, scCube_location[,2])
combined_location$Dataset <- c(rep("SimSpace", nrow(location)), rep("Xenium Reference", nrow(real_location)), rep("scCube", nrow(scCube_location)))
ST_PC <- as.data.frame(t(ST@SpatialPCs))
ST_PC$row <- ST@location[,1]
ST_PC$col <- ST@location[,2]
ST_PC <- left_join(ST_PC, combined_location, by = c("row", "col"))

set.seed(0)
ST_PC <- ST_PC[sample(nrow(ST_PC), nrow(ST_PC)),]

ggplot(data = ST_PC) + 
  geom_point(aes(x = V1, y = V2, color = Cluster), size=1.8, alpha = 0.6, stroke = 0) + 
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33", 
                                "#A65628", "#999999", "#F781BF", "#A65628"),
                     name = "Cell Type") +
  xlab("PC1") +
  ylab("PC2") + 
  theme_bw() + 
  theme(aspect.ratio = 1)

ggplot(data = ST_PC) + 
  geom_point(aes(x = V1, y = V2, color = Dataset), size=1.8, alpha = 0.6, stroke = 0) + 
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33", 
                                  "#A65628", "#999999", "#F781BF", "#A65628"),
                     name = "Source") +
  xlab("PC1") +
  ylab("PC2") + 
  theme_bw() + 
  theme(aspect.ratio = 1)



