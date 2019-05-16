#### CNV-calling benchmark heatmaps ####
#### Christina Bligaard Pedersen ####

#### Set-up ####
# Load packages
library(ComplexHeatmap)
library(RColorBrewer)
library(openxlsx)

# Specify colors from RColorBrewer
lightblue <- brewer.pal(n=12, 'Paired')[1]; darkblue <- brewer.pal(n=12, 'Paired')[2]
lightgreen <- brewer.pal(n=12, 'Paired')[3]; darkgreen <- brewer.pal(n=12, 'Paired')[4]
lightred <- brewer.pal(n=12, 'Paired')[5]; darkred <- brewer.pal(n=12, 'Paired')[6]
lightorange <- brewer.pal(n=12, 'Paired')[7]; darkorange <- brewer.pal(n=12, 'Paired')[8]
lightpurple <- brewer.pal(n=12, 'Paired')[9]; darkpurple <- brewer.pal(n=12, 'Paired')[10]
yellow <- brewer.pal(n=12, 'Paired')[11]; brown <- brewer.pal(n=12, 'Paired')[12]

#out_dir <- "~/Desktop/"
out_dir <- "/home/projects/cu_10047/data/analyzed_samples/project_hierarchy/research_projects/CNV_benchmark/plotting_scripts/"
#in_dir <- "/mnt/computerome_rh/data/analyzed_samples/project_hierarchy/research_projects/CNV_benchmark/concordance/"
in_dir <- "/home/projects/cu_10047/data/analyzed_samples/project_hierarchy/research_projects/CNV_benchmark/concordance/"

# Heatmap function for all true CNVs
make_heatmap <- function(data_frame, color, cl_cols, edge, add_q, title) {
  
  # Get annotation of CNV types and make bar
  if (is.data.frame(add_q)) {
    meta <- qualities[rownames(data_frame),c('Quality_score', 'Manually_checked')]; colnames(meta) <- c('Quality','Man.annot')
    meta_bar <- HeatmapAnnotation(df = meta, col = list(Quality = colorRamp2(c(-1, 0, 1), c(darkblue, "white", darkred)),
                                                        Man.annot = c('Yes' = lightpurple, 'No' = lightorange)), 
                                  gp = gpar(col = "white", lwd = edge), show_annotation_name = T, 
                                  annotation_name_side = 'left', annotation_name_gp = gpar(fontsize = 22),
                                  annotation_legend_param = list(Quality = list(labels_gp = gpar(fontsize = 20),
                                                                                title_gp = gpar(fontsize = 22),
                                                                                border = 'black',
                                                                                grid_height = unit(0.3, "in"),
                                                                                grid_width = unit(0.3, "in")),
                                                                 Man.annot = list(labels_gp = gpar(fontsize = 20),
                                                                                  title_gp = gpar(fontsize = 22),
                                                                                  border = 'black',
                                                                                  grid_height = unit(0.3, "in"),
                                                                                  grid_width = unit(0.3, "in"))))          
  } else {
    meta <- data.frame(gsub('.+_(.+)$', '\\1', rownames(data_frame))); colnames(meta) <- c('Type')
    meta_bar <- HeatmapAnnotation(df = meta, col = list(Type = c('DUP' = lightred, 
                                                                 'DEL' = lightblue,
                                                                 'CNV' = yellow)), 
                                  gp = gpar(col = "white", lwd = edge), show_annotation_name = T, 
                                  annotation_name_side = 'left', annotation_name_gp = gpar(fontsize = 22),
                                  annotation_legend_param = list(Type = list(labels_gp = gpar(fontsize = 20),
                                                                             title_gp = gpar(fontsize = 22),
                                                                             border = 'black',
                                                                             grid_height = unit(0.3, "in"),
                                                                             grid_width = unit(0.3, "in"))))
  }
  
  # Clustering for columns
  if (cl_cols) {
    hc <- hclust(dist(data_frame))
    group <- cutree(hc, k = cl_cols)
    cl <- cluster_within_group(t(data_frame), group)
  } else {
    cl <- T
  }
  
  # Code for actual heatmap
  cnv_heatmap = Heatmap(t(data_frame), column_title_side = 'bottom', column_title = 'True CNVs', 
                        name = "cnv_heatmap", cluster_columns = cl,
                        top_annotation = meta_bar, show_row_dend = T, show_column_dend = T,
                        show_row_names = T, show_column_names = F, column_dend_height = unit(2, 'cm'),
                        row_dend_width = unit(2, 'cm'), row_title_side = 'right',
                        column_title_gp = gpar(fontsize = 22), row_names_gp = gpar(fontsize = 22),
                        col = c("1" = color, "0" = "grey95"), rect_gp = gpar(col = "white", lwd = edge),
                        heatmap_legend_param = list(title = "Match", labels = c('Yes', 'No'), border = 'black', 
                                                    title_gp = gpar(fontsize = 22),
                                                    labels_gp = gpar(fontsize = 20),
                                                    grid_height = unit(0.3, "in"),
                                                    grid_width = unit(0.3, "in")))
  
  # Draw and return plot
  p <- draw(cnv_heatmap, heatmap_legend_side = "right", column_title = paste("CNV calling heatmap at", title), column_title_gp = gpar(fontsize = 32, fontface = "bold"))
  return(cnv_heatmap)
  
}

#### True CNV heatmaps ####

## Process data ----
# Read data
data <- read.table(paste0(in_dir, 'callset.txt'), header = T, sep = '\t')
data$uniq_cnv <- paste0(data$sample, '_', data$cnv)

# Make data into a matrix with tools in columns and individual CNVs in rows
cnv_mat <- data.frame()
for (cnv in unique(data$uniq_cnv)) {
  
  rows <- unique(data[data$uniq_cnv==cnv,])
  
  values <- c()
  for (tool in c(levels(data$tool))) {
    if (tool %in% rows$tool) {
      values <- c(values, as.numeric(!is.na(rows[rows$tool==tool,'value'])))
    } else {
      values <- c(values, NA)
    }
  }
  
  cnv_mat <- rbind(cnv_mat, values)
}
colnames(cnv_mat) <- gsub('GermlineCNVCaller', 'GATK gCNV', levels(data$tool))
rownames(cnv_mat) <- unique(data$uniq_cnv)


# Extract sub-matrices for heatmaps
cnv_mat_wgs <- cnv_mat[grep('GB-WGS-\\d', rownames(cnv_mat)),c('CLC', 'cn.MOPS', 'CNVnator', 'ControlFREEC', 'DELLY', 'GATK gCNV', 'Lumpy', 'Manta')]
cnv_mat_wes <- cnv_mat[grep('GB-WES-\\d', rownames(cnv_mat)),c('CLC', 'cn.MOPS', 'CNVkit', 'CODEX2', 'ExomeDepth', 'GATK gCNV', 'Manta')]
cnv_mat_NA12878_wgs <- cnv_mat[grep('GB-WGS-NA12878', rownames(cnv_mat)),c('CLC', 'cn.MOPS', 'CNVnator', 'ControlFREEC', 'DELLY', 'GATK gCNV', 'Lumpy', 'Manta')]
cnv_mat_NA12878_wes <- cnv_mat[grep('GB-WES-NA12878', rownames(cnv_mat)),c('CLC', 'cn.MOPS', 'CNVkit', 'CODEX2', 'ExomeDepth', 'GATK gCNV', 'Manta')]

# Read quality scores
qualities <- read.xlsx(paste0(in_dir, 'manual_annotation_list.xlsx'), rowNames = T)
qualities$Quality_score <- as.numeric(as.character(qualities$Quality_score))


#### Create heatmaps ####

## WES (8 samples) ----
png(paste0(out_dir, '/CNV_heatmap_WES.png'), width = 20, height = 8, units = 'in', res = 300)
#pdf(paste0(out_dir, '/CNV_heatmap_WES.pdf')', width = 20, height = 8)
make_heatmap(data_frame = cnv_mat_wes, color = darkpurple, cl_cols = F, edge = 1, add_q = qualities, title = "WES level (8 samples)")
dev.off()

## WGS (38 samples) ----
png(paste0(out_dir, '/CNV_heatmap_WGS.png'), width = 20, height = 8, units = 'in', res = 300)
#pdf(paste0(out_dir, '/CNV_heatmap_WGS.pdf'), width = 20, height = 8)
make_heatmap(data_frame = cnv_mat_wgs, color = darkgreen, cl_cols = 20, edge = 0.1, add_q = qualities, title = "WGS level (38 samples)")
dev.off()

## WES (NA12878) ----
png(paste0(out_dir, '/CNV_heatmap_NA12878-WES.png'), width = 8, height = 8, units = 'in', res = 300)
#pdf(paste0(out_dir, '/CNV_heatmap_NA12878-WES.pdf'), width = 8, height = 8)
make_heatmap(data_frame = cnv_mat_NA12878_wes, color = darkred, cl_cols = F, edge = 0.01, add_q = F, title = "WES level (NA12878)")
dev.off()

## WGS (NA12878) ----
png(paste0(out_dir, '/CNV_heatmap_NA12878-WGS.png'), width = 16, height = 8, units = 'in', res = 300)
#pdf(paste0(out_dir, '/CNV_heatmap_NA12878-WGS.pdf'), width = 16, height = 8)
make_heatmap(data_frame = cnv_mat_NA12878_wgs, color = darkred, cl_cols = 20, edge = 0.01, add_q = F, title = "WGS level (NA12878)")
dev.off()


#### Called CNV heatmaps ####
## Process data ----
# Read data
cnv_mat_2 <- read.table(paste0(in_dir, 'all_called.txt'), header = T, sep = '\t')
colnames(cnv_mat_2) <- gsub('GATK_gCNV', 'GATK gCNV', colnames(cnv_mat_2))

cytoscan <- read.table(paste0(in_dir, 'cytoscan_called.txt'), header = T, sep = '\t')
lengths <- as.numeric(gsub('^.+:(\\d+)-(\\d+)_D.+', '\\2', rownames(cnv_mat_2))) - as.numeric(gsub('^.+:(\\d+)-(\\d+)_D.+', '\\1', rownames(cnv_mat_2))) + 2
lengths[is.na(lengths)] <- 2; names(lengths) <- rownames(cnv_mat_2)  # We add two, because the output bed had the start base truncated

# Sanity check
table(rowSums(cnv_mat_2, na.rm = T))

# Split data into WES and WGS sets
cnv_mat_2_wes <- cnv_mat_2[grep('WES', rownames(cnv_mat_2)),c('CLC', 'cn.MOPS', 'CNVkit', 'CODEX2', 'ExomeDepth', 'GATK gCNV', 'Manta')]
cnv_mat_2_wgs <- cnv_mat_2[grep('WGS', rownames(cnv_mat_2)),c('CLC', 'cn.MOPS', 'CNVnator', 'ControlFREEC', 'DELLY', 'GATK gCNV', 'Lumpy', 'Manta')]

#cnv_mat_2_wes <- cnv_mat_2_wes[rowSums(cnv_mat_2_wes, na.rm = T)!=1,]

# Percentages
wgs_percentages <- (table(rowSums(cnv_mat_2_wgs, na.rm = T))/nrow(cnv_mat_2_wgs))*100
wes_percentages <- (table(rowSums(cnv_mat_2_wes, na.rm = T))/nrow(cnv_mat_2_wes))*100


# Heatmap code

# Get annotations
meta <- cbind.data.frame(gsub('.+_(.+)$', '\\1', rownames(cnv_mat_2_wes)), cytoscan[rownames(cnv_mat_2_wes),], lengths[rownames(cnv_mat_2_wes)]<1000); colnames(meta) <- c('Type', 'CytoScan', '>100bp')

meta_bar <- HeatmapAnnotation(df = meta, col = list(Type = c('DUP' = lightred, 
                                                             'DEL' = lightblue),
                                                    CytoScan = c('1' = darkpurple,
                                                                 '0' = lightorange),
                                                    ">100bp" = c('TRUE' = yellow,
                                                                  'FALSE' = darkgreen)), 
                              gp = gpar(col = "white", lwd = 0.01), show_annotation_name = T, 
                              annotation_name_side = 'left', annotation_name_gp = gpar(fontsize = 22),
                              annotation_legend_param = list(Type = list(labels_gp = gpar(fontsize = 20),
                                                                         title_gp = gpar(fontsize = 22),
                                                                         border = 'black',
                                                                         grid_height = unit(0.3, "in"),
                                                                         grid_width = unit(0.3, "in")),
                                                             CytoScan = list(labels_gp = gpar(fontsize = 20),
                                                                             title_gp = gpar(fontsize = 22),
                                                                             border = 'black',
                                                                             grid_height = unit(0.3, "in"),
                                                                             grid_width = unit(0.3, "in"),
                                                                             labels = c('Yes', 'No')),
                                                             ">100bp" = list(labels_gp = gpar(fontsize = 20),
                                                                             title_gp = gpar(fontsize = 22),
                                                                             border = 'black',
                                                                             grid_height = unit(0.3, "in"),
                                                                             grid_width = unit(0.3, "in"),
                                                                             labels = c('Yes', 'No'))))

# Cluster columns
hc <- hclust(dist(cnv_mat_2_wes))
# group <- cutree(hc, k = 20)
# 
# clusters <- cutreeDynamic(hc, distM = as.matrix(dist(x)), method = "tree")
# 
# cl <- as.dendrogram(hc) # Work-around memory error?
# cl <- cluster_within_group(t(cnv_mat_2_wes), group)

# Alternative ordering


# Code for actual heatmap
cnv_heatmap = Heatmap(t(cnv_mat_2_wes), column_title_side = 'bottom', column_title = 'Called CNVs', 
                      name = "cnv_heatmap", cluster_columns = F, column_order = hc$order,
                      top_annotation = meta_bar, show_row_dend = T, show_column_dend = F,
                      show_row_names = T, show_column_names = F,# column_dend_height = unit(2, 'cm'),
                      row_dend_width = unit(2, 'cm'), row_title_side = 'right',
                      column_title_gp = gpar(fontsize = 22), row_names_gp = gpar(fontsize = 22),
                      col = c("1" = darkblue, "0" = "grey95"), rect_gp = gpar(col = "white", lwd = 0.01),
                      heatmap_legend_param = list(title = "Match", labels = c('Yes', 'No'), border = 'black', 
                                                  title_gp = gpar(fontsize = 22),
                                                  labels_gp = gpar(fontsize = 20),
                                                  grid_height = unit(0.3, "in"),
                                                  grid_width = unit(0.3, "in")))

png(paste0(out_dir, '/CNV_heatmap_WES_2.png'), width = 20, height = 8, units = 'in', res = 300)
p <- draw(cnv_heatmap, heatmap_legend_side = "right", column_title = paste("CNV calling heatmap (WES)"), column_title_gp = gpar(fontsize = 32, fontface = "bold"))
dev.off()

