#!/usr/bin/env Rscript

###### Functions for DEA ######

options(warn=-1)

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

createDDS2 <- function(counts, metadata, first.level, ref.level){
  flog.info("########## Create DDS object ###########")

  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = metadata,
                                design = ~ conditions)
  
  
  dds$conditions = factor(dds$conditions, levels = c(first.level, ref.level))
  
  dds <- DESeq(dds, parallel = T)

  return(dds)
}

createRES <- function(dds, first.level, ref.level, pvalue, gtf_file){
  flog.info("########## Create results object ###########")

  res <- results(dds, contrast = c("conditions", first.level, ref.level), alpha = pvalue, parallel = T)
  res_df <- as.data.frame(res)
  res_df <- na.omit(res_df)
  res_df <- res_df[order(res_df$padj, res_df$pvalue, decreasing = F),]
  res_df <- cbind(names = rownames(res_df), res_df)
  res_df$Significance <- ifelse(res_df$padj < pvalue, TRUE, FALSE)
  res_df$gencode =  res_df$names
  tmp = gtf_file[which(gtf_file$gene_id %in% res_df$names),]
  res_df$genes = tmp[match(res_df$names, tmp$gene_id), "gene_name"]
  res_df$genes = make.unique(res_df$genes)
  res_df = na.omit(res_df)
  row.names(res_df) = res_df$genes
  res_df$names = res_df$genes
  res_df$genes <- NULL
  return(res_df)
}

createPCA<- function(rld, first.level, ref.level, condi_col){
  flog.info("########## Create PCA plot ###########")

  rldObject = rld
  pcaData <- plotPCA(rldObject, intgroup="conditions", returnData=TRUE)
  
  percentVar <- round(100 * attr(pcaData, "percentVar"))

  pca_plot = ggplot(pcaData, aes(PC1, PC2, color=conditions, label=name)) +
    scale_color_manual(values = condi_col) +
    geom_point(size=5) +
    ggtitle("PCA plot") +
    theme_bw() +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    guides(color = guide_legend(order = 1), fill = guide_legend(order = 0)) +
    geom_text_repel(aes(label = pcaData$name), size = 6, box.padding = 0.5, max.overlaps = Inf) +
    theme(
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.title = element_text(size = 20, color = "white"),
      legend.key = element_rect(colour = "transparent", fill = "transparent"),
      legend.text = element_text(size = 20, color = "white"),
      axis.text = element_text(angle = 45, hjust = 1, size = 17, color = "white"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "white"),
      axis.title = element_text(size = 23, color = "white"))
  return(pca_plot)
}

createVolcano <- function(res_df, condi_col){
  flog.info("########## Create volcano plot ###########")
  res_df$Significance_reg = ifelse(res_df$Significance, 
                                   ifelse(res_df$log2FoldChange > 0, "Up", "Down"), 
                                   "Not sig.")
  res_df$Significance_reg = factor(res_df$Significance_reg, levels = c("Up", "Down","Not sig."))
  
  res_df$Significance_reg = factor(res_df$Significance_reg, levels = c("Up", "Down","Not sig."))
  
  # Subset significat genes to color in blue
  gen_subset <- subset(res_df, res_df$Significance == TRUE)
  gen_subset <- gen_subset[order(gen_subset$padj),]
  if (dim(gen_subset)[1] > 10){
    gen_subset_short <- gen_subset[1:10,]
  } else {
    gen_subset_short <- gen_subset
  }
  
  res_df$geneLabels = ifelse(res_df$names %in% gen_subset_short$names, TRUE, FALSE)
  color_code = list("Up" = condi_col[1],
                    "Down" = condi_col[2],
                    "Not sig." = "gray")
  
  # Creates volcano plot with the 25 most significant genes labeled
  
  vol = ggplot(data=res_df, aes(x=log2FoldChange, y=-log10(padj), colour=Significance_reg)) +
    #scale_color_manual(values = colorCode) +
    geom_point(size=1.75) +
    xlab("log2 fold change") + ylab("-log10 p-adjusted")+
    ggtitle(paste0("Differential Expression (", names(condi_col)[1], " vs. ", names(condi_col)[2], ")")) + # add conditions to title
    scale_color_manual(values = color_code) +
    theme_bw() +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    #theme(legend.position = "None")+
    geom_text_repel(max.overlaps = Inf, max.time = 1, aes(x=log2FoldChange, y=-log10(padj)),
                    label = ifelse(res_df$geneLabels == TRUE, as.character(res_df$names),""),
                    box.padding = 0.5, show.legend = F, size = 6) +
    theme(
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      # panel.grid.major = element_blank(), # get rid of major grid
      # panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      #legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
      #legend.title = element_text(size = 20, color = "white"),
      legend.title = element_blank(), # remove legend title
      legend.key = element_rect(colour = "transparent", fill = "transparent"),
      legend.text = element_text(size = 20, color = "white"),
      axis.text = element_text(angle = 45, hjust = 1, size = 17, color = "white"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "white"),
      axis.title = element_text(size = 23, color = "white"))
  return(vol)
}

createHeatmap <- function(dds, rld, condi_col, main_color = "RdBu", gtf_file = NA, genes = NA){
  flog.info("########## Create Heatmap of Expression ###########")

  print("Heatmap will be created...")
  colList = list("conditions" = condi_col)
  select <- genes[, "gencode"]
  df <- as.data.frame(colData(dds)["conditions"])
  heat_input = assay(rld)[select,]
  row.names(heat_input) = genes[, "names"]

  ha = HeatmapAnnotation(Condition = df$conditions,
                        col = list(Condition = condi_col),
                        name = "Condition   ",
                        show_annotation_name = F,
                        annotation_height = 2,
                        annotation_width = 2,
                        annotation_legend_param = list(
                          title_gp = gpar(fontsize = 18, col = "white"),
                          labels_gp = gpar(fontsize = 16, col = "white"),
                          title_position = "lefttop-rot"
                        )) # changed from annotation_label = gt_render(c("<span style='color:red'>condition</span>"))
  g = ComplexHeatmap::Heatmap(
    heat_input, 
    name = "Norm. counts",  
    col = hcl.colors(50, main_color),
    cluster_rows = T, 
    cluster_columns = T, 
    show_column_dend = F, 
    top_annotation = ha, 
    show_row_dend = T,
    show_column_names=T,
    column_title = NULL, 
    column_names_rot = 45,
    row_dend_gp = gpar(col = "white"),
    row_names_side = "right", 
    row_names_gp = gpar(fontsize = 15, col = "white"),
    column_names_gp = gpar(fontsize = 20, col = "white"),
    heatmap_legend_param = list(
              title_gp = gpar(fontsize = 18, col = "white"),
              labels_gp = gpar(fontsize = 16, col = "white"),
              legend_height = unit(6, "cm"), 
              grid_width = unit(0.5, "cm"),
              title_position = "lefttop-rot"
              )
  )
  g_draw = draw(g, background = "transparent")

  ha2 = HeatmapAnnotation(Condition = df$conditions,
                         col = list(Condition = condi_col),
                         name = "Condition   ",
                         show_annotation_name = F,
                         annotation_height = 2,
                         annotation_width = 2,
                         annotation_legend_param = list(
                           title_gp = gpar(fontsize = 18, col = "black"),
                           labels_gp = gpar(fontsize = 16, col = "black"),
                           title_position = "lefttop-rot"
                         )) # changed from annotation_label = gt_render(c("<span style='color:red'>condition</span>"))
  
  g2 = ComplexHeatmap::Heatmap(
    heat_input, 
    name = "Norm. counts",  
    col = hcl.colors(50, main_color),
    cluster_rows = T, 
    cluster_columns = T, 
    show_column_dend = F, 
    top_annotation = ha2, 
    show_row_dend = T,
    show_column_names=T,
    column_title = NULL, 
    column_names_rot = 45,
    row_names_side = "right", 
    row_names_gp = gpar(fontsize = 15, col = "black"),
    column_names_gp = gpar(fontsize = 15, col = "black"),
    heatmap_legend_param = list(
              title_gp = gpar(fontsize = 18, col = "black"),
              labels_gp = gpar(fontsize = 16, col = "black"),
              legend_height = unit(6, "cm"), 
              grid_width = unit(0.5, "cm"),
              title_position = "lefttop-rot"
              )
  )

  g_draw2 = draw(g2, background = "transparent")

  return(list("heat" = g_draw, "heat.down" = g_draw2))
}

createSam2Sam <- function(rld){
  flog.info("########## Create Sample to sample distance heatmap ###########")

  sampleDists <- dist(t(assay(rld)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(rld$conditions, rld$Samples, sep="-")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")))(255)

   g = ComplexHeatmap::Heatmap(sampleDistMatrix, 
    name = "Distance", 
                cluster_rows = T, 
                cluster_columns = T, 
                show_column_dend = F,
                show_row_dend = T,
                show_column_names=F, col = colors,
                column_title = NULL, 
                row_names_side = "right", 
                heatmap_legend_param = list(
                  title_gp = gpar(fontsize = 22, col = "white"),
                  labels_gp = gpar(fontsize = 20, col = "white"),
                  legend_direction = "horizontal", 
                  legend_width = unit(8, "cm"), 
                  grid_height = unit(1, "cm"),
                  title_position = "lefttop"),
                row_names_gp = gpar(fontsize = 22, col = "white"),
                row_dend_gp = gpar(col = "white"))
  g_draw = draw(g, background = "transparent", 
                heatmap_legend_side = "bottom", 
                padding = unit(c(2, 2, 2, 30), "mm"))

  g2 = ComplexHeatmap::Heatmap(sampleDistMatrix, 
    name = "Distance", 
                cluster_rows = T, 
                cluster_columns = T, 
                show_column_dend = F,
                show_row_dend = T,
                show_column_names=F, col = colors,
                column_title = NULL, 
                row_names_side = "right", 
                heatmap_legend_param = list(
                  title_gp = gpar(fontsize = 12, col = "black"),
                  labels_gp = gpar(fontsize = 10, col = "black"),
                  legend_direction = "horizontal",
                  legend_width = unit(4, "cm"),
                  grid_height = unit(0.5, "cm"),
                  title_position = "lefttop"),
                row_names_gp = gpar(fontsize = 12, col = "black"),
                row_dend_gp = gpar(col = "black"))
  g_draw2 = draw(g2, background = "transparent", 
                 heatmap_legend_side = "bottom", 
                 padding = unit(c(2, 2, 2, 2), "mm"))

  return(list(g_draw, g_draw2))
}

run_preprocessing_dea <- function(meta.file, counts.file, condition.col, first.level, ref.level, pvalue, gtf_file){
  flog.info("########## Differential Expression Analysis ###########")
  counts = counts.file
  
  metadata = meta.file
  row.names(metadata) <- metadata$Samples
  
  missingSampleInfos = colnames(counts)[-which(colnames(counts) %in% metadata$Samples)]
  if (length(missingSampleInfos) > 0){
    print(paste0("No metadata found for the following samples: ", paste(missingSampleInfos, collapse = ",")))
    print(paste0(">>>>> Counts will be excluded!"))
    
  } else {
    print("All required information is included!")
  }
  metadata = metadata[which(row.names(metadata) %in% colnames(counts)),]
  
  colnames(metadata)[colnames(metadata) == condition.col] = "conditions"
  
  if (is.na(first.level) & is.na(ref.level)){
    first.level = unique(metadata$conditions)[1]
    ref.level = unique(metadata$conditions)[2]
    
  } else if (is.na(first.level)){
    first.level = unique(metadata$conditions[!(metadata$conditions == ref.level)])[1]
  } else {
    ref.level = unique(metadata$conditions[!(metadata$conditions == first.level)])[1]
  }
  print(metadata)
  print(c(first.level, ref.level))
  metadata = metadata[which(metadata$conditions %in% c(first.level, ref.level)),]
  metadata = metadata[which(row.names(metadata) %in% intersect(row.names(metadata), colnames(counts))),]
  
  counts = counts[, match(row.names(metadata), colnames(counts))]
  
  dds = createDDS2(counts, metadata, first.level, ref.level)
  rld = rlog(dds)
  res_df = createRES(dds, first.level, ref.level, pvalue, gtf_file)
  print(head(res_df))
  
  deaProcess = list(res_df = res_df, 
                    rld = rld, 
                    dds = dds, 
                    counts = counts, 
                    metadata = metadata, 
                    first.level = first.level, 
                    ref.level = ref.level)
  
}

save_rds <- function(deaResults, output.dir){
  
  flog.info("Saving results to the rds object: deaResults.RDS")
  print(names(deaResults))
  print(paste0(output.dir, "deaResults.rds"))
  saveRDS(deaResults, paste0(output.dir, "deaResults.rds"))
  
}
