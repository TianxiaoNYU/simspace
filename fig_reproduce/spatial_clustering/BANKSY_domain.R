library(Banksy)

library(SummarizedExperiment)
library(SpatialExperiment)
library(scuttle)

library(scater)
library(cowplot)
library(ggplot2)
library(dplyr)

# Load the data
{
  rawcount <- read.csv('../../fig_data/clustering/simspace_count.csv', 
                       header = T,
                       row.names = 1)
  location <- read.csv('../../fig_data/clustering/simspace_metadata.csv',
                       header = T, )
  rawcount <- t(rawcount)
  # location <- as.matrix(location)
  
  rownames(location) <- colnames(rawcount)
  simspace_domain <- location
  new_location <- location[,c("row", "col", 'fitted_celltype')]
  # new_location <- as.data.frame(new_location)
  
  ## scCube data loading
  {
    scCube_rawcount <- read.csv('../../fig_data/clustering/scCube_count.csv', 
                         header = T,
                         row.names = 1)
    scCube_location <- read.csv('../../fig_data/clustering/scCube_metadata.csv',
                         header = T)
    scCube_location$point_x <- (scCube_location$point_x  - min(scCube_location$point_x )) / (max(scCube_location$point_x ) - min(scCube_location$point_x )) * 100
    scCube_location$point_y <- (scCube_location$point_y  - min(scCube_location$point_y )) / (max(scCube_location$point_y ) - min(scCube_location$point_y )) * 100

    scCube_location$point_x <- scCube_location$point_x + 100
    scCube_location <- as.matrix(scCube_location)
    
    rownames(scCube_location) <- colnames(scCube_rawcount)
    new_scCube_location <- scCube_location[,c("point_x", "point_y", 'Cell_type')]
    new_scCube_location <- as.data.frame(new_scCube_location)
    colnames(new_scCube_location) <- c("row", "col", 'state_rank')
  }
  
  real_rawcount <- read.csv('../../data/reference_count.csv',
                            header = T,
                            row.names = 1)
  real_meta <- read.csv('../../data/reference_metadata.csv',
                        header = T)
  cell_count <- sort(table(real_meta$Cluster), decreasing = T)
  real_meta$state_rank <- 0
  for (i in 1:length(cell_count)) {
    real_meta$state_rank[real_meta$Cluster == names(cell_count)[i]] <- i 
  }
  real_location <- as.data.frame(real_meta[,c("x_centroid", "y_centroid", "Cluster")])
  # normalize the coordinates in real_location so that they sit between 0 and 100 for each column seperately
  real_location$x_centroid <- (real_location[,1] - min(real_location[,1])) / (max(real_location[,1]) - min(real_location[,1])) * 100
  real_location$y_centroid <- (real_location[,2] - min(real_location[,2])) / (max(real_location[,2]) - min(real_location[,2])) * 100
  colnames(real_location) <- c("row", "col")
  real_location$row <- real_location$row + 200
  rownames(real_location) <- paste0("real_", seq(1, nrow(real_location)))
  colnames(real_rawcount) <- rownames(real_location)
  real_rawcount <- as.matrix(real_rawcount)
  
  
  
  combined_rawcount <- cbind(rawcount, scCube_rawcount, real_rawcount)
  new_location <- as.matrix(new_location)
  real_location <- as.matrix(real_location)
  new_scCube_location <- as.matrix(new_scCube_location)
  combined_location <- rbind(new_location, new_scCube_location, real_location)
  combined_location <- as.data.frame(combined_location)
  combined_location$Dataset <- c(rep("SimSpace", nrow(location)), rep("scCube", nrow(scCube_location)), rep("Xenium Reference", nrow(real_location)))
  
  combined_rawcount <- as.matrix(combined_rawcount)
  combined_location <- as.matrix(combined_location)
}

combined_location <- as.data.frame(combined_location)
combined_location$row <- as.numeric(combined_location$row)
combined_location$col <- as.numeric(combined_location$col)

se <- SpatialExperiment(assay = list(counts = combined_rawcount), 
                        spatialCoords = as.matrix(combined_location[,c(1,2)]))

# Normalization to mean library size
se <- computeLibraryFactors(se)
aname <- "normcounts"
assay(se, aname) <- normalizeCounts(se, log = FALSE)

lambda <- c(0, 0.2, 0.5, 0.8)
k_geom <- c(15, 30)

se <- Banksy::computeBanksy(se, assay_name = aname, compute_agf = TRUE, k_geom = k_geom)

set.seed(1)
se <- Banksy::runBanksyPCA(se, use_agf = TRUE, lambda = lambda)
se <- Banksy::runBanksyUMAP(se, use_agf = TRUE, lambda = lambda)
se <- Banksy::clusterBanksy(se, use_agf = TRUE, lambda = lambda, resolution = 1.2)
se <- Banksy::connectClusters(se)

cnames <- colnames(colData(se))
cnames <- cnames[grep("^clust", cnames)]
colData(se) <- cbind(colData(se), spatialCoords(se))

tmp <- as.data.frame(se@int_colData)
combined_location$PCA1 <- tmp$reducedDims.PCA_M1_lam0.5.PC1
combined_location$PCA2 <- tmp$reducedDims.PCA_M1_lam0.5.PC2
combined_location$UMAP1 <- tmp$reducedDims.UMAP_M1_lam0.5.V1
combined_location$UMAP2 <- tmp$reducedDims.UMAP_M1_lam0.5.V2
combined_location$UMAP1_cluster <- se@colData$clust_M1_lam0.5_k50_res1.2

combined_location$domain_cluster <- "Domain 1"
combined_location$domain_cluster[combined_location$UMAP1_cluster %in% c("6")] <- "Domain 2"
combined_location$domain_cluster[combined_location$UMAP1_cluster %in% c("7")] <- "Domain 3"
combined_location$domain_cluster[combined_location$UMAP1_cluster %in% c("1", "5", "8")] <- "Domain 4"

combined_location$Tumor_label <- "Tumor"
combined_location$Tumor_label[combined_location$UMAP1_cluster %in% c("6")] <- "Stromal"

# set.seed(0)
# combined_location <- combined_location[sample(nrow(combined_location), nrow(combined_location)),]

ggplot(data = combined_location) + 
  geom_point(aes(x = row, y = col, color = UMAP1_cluster, shape = Dataset), size=1.7) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33", 
                                "#A65628", "#999999", "#F781BF", "#C65628", "#008695FF", "#4B4B8FFF"),
                     name = "Cell Cluster") +
  scale_shape_manual(values = c(16, 17, 15), name = "Source") +
  xlab("X") + 
  ylab("Y") +
  theme_bw() + 
  theme(aspect.ratio = 0.33) 

ggplot(data = combined_location) + 
  geom_point(aes(x = UMAP1, y = UMAP2, color = UMAP1_cluster), size=1.8, alpha = 0.6, stroke = 0) + 
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33", 
                                "#A65628", "#999999", "#F781BF", "#C65628", "#008695FF", "#4B4B8FFF"),
                     name = "Cell Cluster") +
  theme_bw() + 
  theme(aspect.ratio = 1)

ggplot(data = combined_location) + 
  geom_point(aes(x = UMAP1, y = UMAP2, color = domain_cluster), size=1.8, alpha = 0.75, stroke = 0) + 
  theme_bw() + 
  theme(aspect.ratio = 1)

ggplot(data = combined_location) + 
  geom_point(aes(x = row, y = col, color = domain_cluster, shape = Dataset), size=1.8, alpha = 0.75, stroke = 0) + 
  scale_shape_manual(values = c(16, 17, 15), name = "Source") +
  theme_bw() + 
  xlab("X") + 
  ylab("Y") +
  theme(aspect.ratio = 0.33) 

#### Domain analysis
bank_simspace <- combined_location %>%
  filter(Dataset == "SimSpace")
simspace_domain <- simspace_domain[rownames(bank_simspace), ]
bank_simspace$niche <- simspace_domain$niche
table(bank_simspace$niche, bank_simspace$domain_cluster)
# write.table(bank_simspace,
#             '/Users/zhaotianxiao/Library/CloudStorage/Dropbox/FenyoLab/Project/Spatialsim/SimSpace/fig_data/clustering/BANKSY_ss_domain.csv',
#             sep = ",", quote = F)
  
ggplot(data = bank_simspace) + 
  geom_point(aes(x = row, y = col, color = domain_cluster), size=1.8, alpha = 0.75, stroke = 0) + 
  theme_bw() + 
  xlab("X") + 
  ylab("Y") +
  theme(aspect.ratio = 1) 
ggplot(data = bank_simspace) + 
  geom_point(aes(x = row, y = col, color = UMAP1_cluster), size=1.8, alpha = 0.75, stroke = 0) + 
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33", 
                                "#A65628", "#999999", "#F781BF", "#C65628", "#008695FF", "#4B4B8FFF"),
                     name = "Cell Cluster") +
  theme_bw() + 
  xlab("X") + 
  ylab("Y") +
  theme(aspect.ratio = 1) 



bank_xenium <- combined_location %>%
  filter(Dataset == "Xenium Reference")
ggplot(data = bank_xenium) + 
  geom_point(aes(x = row, y = col, color = domain_cluster), size=1.8, alpha = 0.75, stroke = 0) + 
  theme_bw() + 
  xlab("X") + 
  ylab("Y") +
  theme(aspect.ratio = 1) 
ggplot(data = bank_xenium) + 
  geom_point(aes(x = UMAP1, y = UMAP2, color = domain_cluster), size=1.8, alpha = 0.75, stroke = 0) + 
  theme_bw() + 
  xlab("X") + 
  ylab("Y") +
  theme(aspect.ratio = 1) 
# write.table(bank_xenium,
#             '/Users/zhaotianxiao/Library/CloudStorage/Dropbox/FenyoLab/Project/Spatialsim/SimSpace/fig_data/clustering/BANKSY_xenium_domain.csv',
#             sep = ",", quote = F)
