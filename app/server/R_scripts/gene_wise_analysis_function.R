#!/usr/bin/env Rscript

###### Gene-wise analyses  #####
options(warn=-1)

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

getGeneSymbolFromGTF <- function(gtf.file, output.dir){
  gtf.gr = rtracklayer::import(gtf.file) # creates a GRanges object
  gtf.df = as.data.frame(gtf.gr)
  genes = unique(gtf.df[ ,c("gene_id","gene_name", "transcript_id", "transcript_name")])
  return(genes)
}
createDDS <- function(counts.file, meta.file, condition.col, first.level, ref.level){
  print("########## Normalization of counts plots ###########")
  counts = counts.file

  metadata = meta.file
  row.names(metadata) <- metadata$Samples
  
  missingSampleInfos = colnames(counts)[-which(colnames(counts) %in% metadata$Samples)]
  if (length(missingSampleInfos) > 0){
    print(paste0("No metadata found for the following samples: ", paste(missingSampleInfos, collapse = ",")))
    print(paste0(">>>>> Counts will be excluded!"))
    
  } else {
    print("All required informations are included!")
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
  
  metadata = metadata[which(metadata$conditions %in% c(first.level, ref.level)),]
  metadata = metadata[which(row.names(metadata) %in% intersect(row.names(metadata), colnames(counts))),]
  counts = counts[, match(row.names(metadata), colnames(counts))]
  
  x = all(row.names(metadata) == colnames(counts))
  if (x){
    print("Starting normalization ...")
  } else {
    print("Sample names do not match between counts and metadata!!")
  }
  
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = metadata,
                                design = ~ conditions)
  
  
  dds$conditions = factor(dds$conditions, levels = c(first.level, ref.level))
  
  dds <- DESeq(dds, parallel = T)
  norm_counts <- counts(dds, normalize = TRUE)
  
  return(list(counts, norm_counts, metadata))
}
createCountsPlot <- function(normCounts, genes, metaTab, genes.tab, gtitle, outName, outDir = ".", condi_cols, download = F, ylabel = "Counts"){
  print("########## Create counts plots ###########")
  
  normmutGenes <- normCounts
  print(genes)
  mutGenes <- as.data.frame(normmutGenes[as.character(genes[,1]),])
  mutGenes$genes <- rownames(mutGenes)
  meltedmutGenes <- melt(mutGenes)
  colnames(meltedmutGenes) <- c("gene", "samplename", "normalized_counts")
  meltedmutGenes <- merge(meltedmutGenes, metaTab, by.x = "samplename", by.y = "row.names")
 
  meltedmutGenes_all = merge.data.frame(meltedmutGenes, genes, by.x = "gene", by.y = "gene_id", all.x = T)
  if (!download){
    x1 = ggplot(meltedmutGenes_all) +
      geom_boxplot(aes(x = gene_name, y = normalized_counts, fill = conditions), color = "white") +
      #geom_point(aes(x = gene_name, y = normalized_counts, color = conditions), position = position_jitter(w=0.1, h=0), size=4)+
      scale_y_log10(oob = scales::squish_infinite) +
      scale_fill_manual(values = condi_cols) +
      
      ylab(ylabel) +
      xlab("Genes") +
      ggtitle(gtitle) +
      theme_bw() +
      theme(
        panel.background = element_rect(fill = "transparent", color = "white"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "white"), # get rid of major grid
        panel.grid.minor = element_line(size = 0.2, linetype = 'solid', colour = "white"), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
        legend.title = element_text(size = 20, color = "white"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.text = element_text(size = 20, color = "white"),
        axis.text = element_text(angle = 45, hjust = 1, size = 17, color = "white"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "white"),
        axis.title = element_text(size = 23, color = "white"))
    
    x2 <- ggplot(meltedmutGenes_all) +
      #geom_boxplot(aes(x = gene_name, y = normalized_counts, fill = conditions)) +
      geom_point(aes(x = gene_name, y = normalized_counts, color = conditions), position = position_jitter(w=0.1, h=0), size=4)+
      scale_y_log10(oob = scales::squish_infinite) +
      scale_color_manual(values = condi_cols) +
      #scale_fill_manual(values = safe_colorblind_palette[c(3,11)]) +
      
      ylab(ylabel) +
      xlab("Genes") +
      ggtitle(gtitle) +
      theme_bw() +
      theme(
        panel.background = element_rect(fill = "transparent", color = "white"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "white"), # get rid of major grid
        panel.grid.minor = element_line(size = 0.2, linetype = 'solid', colour = "white"), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
        legend.title = element_text(size = 20, color = "white"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.text = element_text(size = 20, color = "white"),
        axis.text = element_text(angle = 45, hjust = 1, size = 17, color = "white"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "white"),
        axis.title = element_text(size = 23, color = "white"))
  } else {
       x1 = ggplot(meltedmutGenes_all) +
      geom_boxplot(aes(x = gene_name, y = normalized_counts, fill = conditions), color = "black") +
      scale_y_log10(oob = scales::squish_infinite) +
      scale_fill_manual(values = condi_cols) +
      
      ylab(ylabel) +
      xlab("Genes") +
      ggtitle(gtitle) +
      theme_bw() +
      theme(
        legend.title = element_text(size = 20, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        axis.text = element_text(angle = 45, hjust = 1, size = 17, color = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "black"),
        axis.title = element_text(size = 23, color = "black"))
    
    x2 <- ggplot(meltedmutGenes_all) +
      geom_point(aes(x = gene_name, y = normalized_counts, color = conditions), position = position_jitter(w=0.1, h=0), size=4)+
      scale_y_log10(oob = scales::squish_infinite) +
      scale_color_manual(values = condi_cols) +
      ylab(ylabel) +
      xlab("Genes") +
      ggtitle(gtitle) +
      theme_bw() +
      theme(
        legend.title = element_text(size = 20, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        axis.text = element_text(angle = 45, hjust = 1, size = 17, color = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "black"),
        axis.title = element_text(size = 23, color = "black"))

  }
  return(list("Points" = x2, "Boxplot" = x1))
}
TEA <- function(counts, norm_counts, genes.list, metadata, pvalue, output.dir, condi_cols){
  print("###################################")
  print("##         Starting TEA          ##")
  print("###################################")
  norm_counts_plot = createCountsPlot(
    normCounts = norm_counts, 
    genes = genes.list, 
    metaTab = metadata, 
    gtitle = "Normalized counts", 
    outName = "normalized", 
    outDir = output.dir, condi_cols = condi_cols, 
    ylabel = "Normalized read counts")

  norm_counts_plot.download = createCountsPlot(
    normCounts = norm_counts, 
    genes = genes.list, 
    metaTab = metadata, 
    gtitle = "Normalized counts", 
    outName = "normalized", 
    outDir = output.dir, condi_cols = condi_cols, download = T,
    ylabel = "Normalized read counts")
  
  counts_plot = createCountsPlot(
    normCounts = counts, 
    genes = genes.list, 
    metaTab = metadata, 
    gtitle = "Raw counts", 
    outName = "raw", 
    outDir = output.dir, condi_cols = condi_cols,
    ylabel = "Raw read counts")

  counts_plot.download = createCountsPlot(
    normCounts = counts, 
    genes = genes.list, 
    metaTab = metadata, 
    gtitle = "Raw counts", 
    outName = "raw", 
    outDir = output.dir, condi_cols = condi_cols, download = T,
    ylabel = "Raw read counts")


  p1 = ggarrange(plotlist = list(counts_plot[["Points"]], norm_counts_plot[["Points"]]), nrow = 1, ncol = 2, common.legend = TRUE)
  p2 = ggarrange(plotlist = list(counts_plot[["Boxplot"]], norm_counts_plot[["Boxplot"]]), nrow = 1, ncol = 2, common.legend = TRUE)
  
  p1.download = ggarrange(plotlist = list(counts_plot.download[["Points"]], norm_counts_plot.download[["Points"]]), nrow = 1, ncol = 2, common.legend = TRUE)
  p2.download = ggarrange(plotlist = list(counts_plot.download[["Boxplot"]], norm_counts_plot.download[["Boxplot"]]), nrow = 1, ncol = 2, common.legend = TRUE)

  return(list("Dotplot" = p1, "Boxplot" = p2, "Dotplot.down" = p1.download, "Boxplot.down" = p2.download))
}
geneBodyCov.plot <- function(gB_results, geneOfInterest_in, metadata, condi_cols){
  
  theme_set(theme_light())

  geneOfInterest_ID = geneOfInterest_in$gene_id
  geneOfInterest_name = geneOfInterest_in$gene_name
  geneOfInterest = paste0(geneOfInterest_name, "\n(", geneOfInterest_ID, ")")
  geneBodyCov = gB_results %>%
    gather(key = "Position", "Value", -Percentile)  %>%
    group_by(Percentile) %>%
    mutate(Position = as.numeric(gsub("X", "", Position)), PercVal = ifelse(Value == 0, 0, Value/sum(Value)))
  n_samples = length(unique(geneBodyCov$Percentile))
  g = ggplot(geneBodyCov, aes(x = Position, y = PercVal, color = Percentile)) +
    geom_smooth(se = FALSE) +
    ylim(0,max(geneBodyCov$PercVal)) +
    ggtitle(geneOfInterest) +
    scale_color_manual("Samples", values = safe_colorblind_palette[c(1:n_samples)]) +
    scale_x_continuous(breaks = c(0, 50, 100), labels = c("3'", "mid gene", "5'")) + # change position labels from 0 to 100 to 3' to 5'
    ylab("Relative Coverage (%)") +
    theme(
    panel.grid.minor.x = element_blank(), # remove minor grid lines from plot
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    # panel.grid.major = element_blank(), # get rid of major grid
    # panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    #legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.title = element_text(size = 20, color = "white"),
    legend.key = element_rect(colour = "transparent", fill = "transparent"),
    legend.text = element_text(size = 20, color = "white"),
    axis.text.y = element_text(angle = 45, hjust = 1, size = 17, color = "white"),
    axis.text.x = element_text(size = 17, color = "white"), # removed: angle = 45, hjust = 1, due to the new discrete axis labels
    plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "white"),
    axis.title = element_text(size = 23, color = "white"),
    )

  geneBodyCov = geneBodyCov %>% 
    left_join(metadata, by = c("Percentile" = "Samples"))
  
 g_cond = ggplot(geneBodyCov, aes(x = Position, y = PercVal, color = Condition)) +
    geom_smooth(se = FALSE) +
    ylim(0,max(geneBodyCov$PercVal)) +
    ggtitle(geneOfInterest) +
    scale_color_manual("Condition", values = condi_cols) +
    scale_x_continuous(breaks = c(0, 50, 100), labels = c("3'", "mid gene", "5'")) + # change position labels from 0 to 100 to 3' to 5'
    ylab("Relative Coverage (%)") +
    theme(
      panel.grid.minor.x = element_blank(), # remove minor grid lines from plot
      panel.background = element_rect(fill = "transparent", color = "white"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                      colour = "white"), # get rid of major grid
      panel.grid.minor = element_line(size = 0.2, linetype = 'solid', colour = "white"), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
      legend.title = element_text(size = 20, color = "white"),
      legend.key = element_rect(colour = "transparent", fill = "transparent"),
      legend.text = element_text(size = 20, color = "white"),
      axis.text.y = element_text(angle = 45, hjust = 1, size = 17, color = "white"),
      axis.text.x = element_text(size = 17, color = "white"),  # removed: angle = 45, hjust = 1, due to the new discrete axis labels
      plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "white"),
      axis.title = element_text(size = 23, color = "white"))



  return(list("samples" = g, "condition" = g_cond))

}