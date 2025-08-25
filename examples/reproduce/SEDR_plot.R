library(ggplot2)

combined_location <- read.table(
  '../../fig_data/clustering/SEDR_umap_embedding.csv',
  header = TRUE, 
  sep = ',',
  row.names = 1)
set.seed(0)
combined_location <- combined_location[sample(nrow(combined_location), nrow(combined_location)),]


ggplot(data = combined_location) + 
  geom_point(aes(x = UMAP1, y = UMAP2, color = Cluster), size=1.8, alpha = 0.6, stroke = 0) + 
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33", 
                                "#A65628", "#999999", "#F781BF", "#A65628"),
                     name = "Cell Type") +
  theme_bw() + 
  theme(aspect.ratio = 1)

ggplot(data = combined_location) + 
  geom_point(aes(x = UMAP1, y = UMAP2, color = Dataset), size=1.8, alpha = 0.6, stroke = 0) + 
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"),
                     name = "Source") +
  theme_bw() + 
  theme(aspect.ratio = 1)
