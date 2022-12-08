

DRIM_seq_prep <- function(table = count.table, run.dir = csv.dir, samps = metadata, condition_col = "Condition", first.level = "a2d3-OE", ref.level = "Ctrl", gtf_file = gtf,cores = 4){
      
    ###############################################################################################
    #                                                                                             #
    #                                                                                             #
    #                                                                                             #
    #                             Read count file, gtf and metadata                               #
    #                                                                                             #
    #                                                                                             #
    #                                                                                             #
    #                                                                                             #
    ###############################################################################################
    
    txdb.filename <- str_split(gtf_file, "/")
    txdb.filename <- as.vector(txdb.filename[[1]])[length(txdb.filename[[1]])]
    txdb.filename <- paste0(run.dir,txdb.filename)
    samps["sample_id"] = samps$Samples
    samps["condition"] = samps[condition_col]
    #print(samps$condition)
    
    samps <- samps[which(samps$condition %in% c(first.level,ref.level)),]
    print(samps)
      
    head(table)
    table = table[,which(colnames(table) %in% samps$Samples)]
    #head(table)
    #txdb <- makeTxDbFromGFF(gtf)
    #print(txdb)
    #saveDb(txdb, txdb.filename)
    #txdb <- loadDb(txdb.filename)
    
    out <- tryCatch(
      {txdb <- loadDb(txdb.filename)
      x = T
      },
      error = function(e){
        x = F
      },
      finally = {
      })
    
    
    #print(out)
    
    if (out){
      txdb <- loadDb(txdb.filename)
    }
    else {
      txdb <- makeTxDbFromGFF(gtf_file)
      saveDb(txdb, txdb.filename)
    }
    
    txdf <- AnnotationDbi::select(txdb, keys(txdb, "GENEID"),"TXNAME", "GENEID")
    #print(txdf)
    tab <- table(txdf$GENEID)
    txdf$ntx <- tab[match(txdf$GENEID, names(tab))]
    #print(txdf$ntx)
    
    
    all(rownames(table) %in% txdf$TXNAME)
    txdf <- txdf[match(rownames(table),txdf$TXNAME),]
    all(rownames(table) == txdf$TXNAME)
    
    counts <- data.frame(gene_id=txdf$GENEID,
                         feature_id=txdf$TXNAME,
                         table)
    #counts <- counts[which(counts$gene_id == goi_id),]
    #print("Counts")
    #head(counts)
    
  
  param = BiocParallel::SerialParam()
  d <- dmDSdata(counts=counts, samples=samps)
  n <- length(samps$sample_id)
  d <- dmFilter(d,
                min_samps_feature_expr=as.integer(n/2), min_feature_expr=10,
                #              min_samps_feature_prop=int(n.small/1.5), min_feature_prop=0.1,
                min_samps_gene_expr=(n/2), min_gene_expr=50)
  
  table(table(counts(d)$gene_id))
  design_full <- model.matrix(~condition, data=DRIMSeq::samples(d))
  #print(design_full)
  set.seed(1)
  system.time({
    d <- dmPrecision(d, design=design_full, BPPARAM = param)
    d <- dmFit(d, design=design_full, BPPARAM = param)
    d <- dmTest(d, coef=colnames(design_full)[2], BPPARAM = param)
  })
  print("DRIM SEQ MADE IT")
  #print(counts(d))
  list = list()
  list$counts = counts
  list$drim = d
  list$samps = samps
  list$txdf = txdf
  return(list)
}


###############################################################################################
#                                                                                             #
#                                                                                             #
#                                                                                             #
#                                             DTU                                             #
#                                                                                             #
#                                                                                             #
#                                                                                             #
#                                                                                             #
###############################################################################################


DTU_special <- function(d_list = tryout, condition_col = "Condition", first.level = "a2d3-OE", ref.level = "Ctrl", goi_id = "ENSG00000111640.15",gtf_tab = gtf_table, cores = 4){
    
    d = d_list$drim
    counts = d_list$counts
    samps = d_list$samps
    
    res <- DRIMSeq::results(d)
    head(res)
    
    res.txp <- DRIMSeq::results(d, level="feature")
    head(res.txp)
    
    no.na <- function(x) ifelse(is.na(x), 1, x)
    res$pvalue <- no.na(res$pvalue)
    res.txp$pvalue <- no.na(res.txp$pvalue)
    
    pScreen <- res$pvalue
    strp <- function(x) substr(x,1,15)
    names(pScreen) <- strp(res$gene_id)
    pConfirmation <- matrix(res.txp$pvalue, ncol=1)
    rownames(pConfirmation) <- strp(res.txp$feature_id)
    
    tx2gene <- res.txp[,c("feature_id", "gene_id")]
    for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])
    goi_df = res.txp[which(res.txp$gene_id == goi_id),]
    goi_df_merged = merge(goi_df,counts(d), by = "feature_id")
    idx <- which(res$gene_id == goi_id)
    #plotProportions(d, res$gene_id[idx], "condition")
    selected_d = counts(d)[which(counts(d)$gene_id == res$gene_id[idx]),]
    samples = samps$sample_id
    feature_ids = unique(selected_d$feature_id)
    plot_dataframe_d = data.frame()
    sum_df = data.frame(sample = as.character(), sum = as.numeric())
    for (i in samples){
      sum = sum(selected_d[i])
      temp_sum = cbind(i,sum)
      colnames(temp_sum) = c("sample", "sum")
      sum_df = rbind(sum_df,temp_sum)
    }
    print(sum_df)
    for (i in samples){
      for (j in feature_ids){
        counts_df = as.numeric(selected_d[which(selected_d$feature_id == j),i])
        feature_id = j
        sample_name = i
        sum_column = as.numeric(sum_df[which(sum_df$sample == i),"sum"])
        percentage = as.numeric(counts_df / sum_column)
        Condition = as.character(samps[which(samps$Samples == i), "Condition"])
        temp_df = cbind(feature_id,counts_df, sample_name, Condition, sum_column, percentage)
        plot_dataframe_d = rbind(plot_dataframe_d,temp_df)
      }
      
      
    }
    plot_dataframe_d$counts_df = as.numeric(plot_dataframe_d$counts_df)
    plot_dataframe_d$percentage = as.numeric(plot_dataframe_d$percentage)
    plot_dataframe_d$Significance = res.txp[which(res.txp$gene_id == goi_id),]$adj_pvalue < 0.05
    plot_dataframe_d$Significance[is.na(plot_dataframe_d$Significance)] <- FALSE
    goi_name = unique(gtf_tab$gene_name[gtf_tab$gene_id == goi_id])
    print(goi_name)
    #plot_dataframe_d  
  
    if(dim(plot_dataframe_d)[1] > 0){
    bp <- ggboxplot(plot_dataframe_d, "feature_id", "percentage",
                  color = "Condition", fill = "Significance", add = "jitter", add.params = list(size = 3, alpha = 1)) +
    color_palette(palette = "jco")+
    fill_palette(palette =  c("steelblue4","indianred3"))+
    xlab("Feature ID") +
    ylab("Percentage") +
    ggtitle(goi_name) +
    theme(
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      #legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
      legend.title = element_text(size = 20, color = "white"),
      legend.key = element_rect(colour = "transparent", fill = "transparent"),
      legend.text = element_text("Condition",size = 8, color = "white"),
      axis.line = element_line(color = "white"),
      axis.text = element_text(angle = 45, hjust = 1, size = 10, color = "white"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "white"),
      axis.title = element_text(size = 14, color = "white"))
  
  bp <- ggpar(bp,legend = "right",
              font.legend = c(14))
  output_list = list()
  output_list$bp = bp
  return(output_list)
    }
  else{
    output_list = list()
    x = NA
    output_list$bp = x
    return(output_list)
  }
}

  
















###############################################################################################
#                                                                                             #
#                                                                                             #
#                                                                                             #
#                                           DTE                                               #
#                                                                                             #
#                                                                                             #
#                                                                                             #
#                                                                                             #
###############################################################################################


DTE_general <- function(d_list, condition_col = "Condition", first.level = "Hct116", ref.level = "MCF7", samps = metadata, gtf_table, cores = 4){
  ###############################################################################################
  #                                                                                             #
  #                                                                                             #
  #                                                                                             #
  #                                   DTE with DEXSeq                                           #
  #                                                                                             #
  #                                                                                             #
  #                                                                                             #
  #                                                                                             #
  ###############################################################################################
  output_list = list()
  d = d_list$drim
  counts = d_list$counts
  samps= d_list$samps

  sample.data <- DRIMSeq::samples(d)
  print(sample.data)
  count.data <- round(as.matrix(counts(d)[,-c(1:2)]))
  print(count.data)
  dxd <- DEXSeqDataSet(countData=count.data,
                       sampleData=sample.data,
                       design=~sample + exon + condition:exon,
                       featureID=counts(d)$feature_id,
                       groupID=counts(d)$gene_id)
  #print(dxd)
  
  system.time({
    dxd <- estimateSizeFactors(dxd)
    dxd <- estimateDispersions(dxd, quiet=TRUE)
    dxd <- testForDEU(dxd, reducedModel=~sample + exon)
    dxd <- estimateExonFoldChanges( dxd, fitExpToVar="condition")
  })
  dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)
  output_list = list()
  dxr_df = as.data.frame(na.omit(dxr))
  print(colnames(dxr_df))
  output_list$dxr = dxr_df
  qval <- perGeneQValue(dxr)
  dxr.g <- data.frame(gene=names(qval),qval)
  
  columns <- c("featureID","groupID","pvalue")
  dxr <- as.data.frame(dxr[,columns])
  #print(head(dxr))
  dxr_df = dxr_df[order(dxr_df$padj, decreasing = FALSE), ]
  
  counter=0
  data_to_label = c()
  for (i in dxr_df$padj){
    if (counter < 10){
      data_to_label = c(data_to_label,TRUE)
    }
    else{
      data_to_label = c(data_to_label,FALSE)
    }
    counter = counter + 1
  }
  print("Line X passed")
  
  dxr_df$label = data_to_label
  dxr_df

  regex = c('log2fold') 
  column_list = as.vector(colnames(dxr_df))
  log_2_fold_name = as.character(column_list[which(grepl(regex,column_list))])
  dxr_df["log_2_fold_Change"] = dxr_df[log_2_fold_name]

  
  print("Line Y passed")
  
  dxr_df$Significance = dxr_df$padj < 0.05
  dxr_df$Significance_reg = ifelse(dxr_df$Significance, 
                                   ifelse(dxr_df$log_2_fold_Change > 0, "Up", "Down"), 
                                   "Not sig.")
  
  dxr_df$Significance_reg = factor(dxr_df$Significance_reg, levels = c("Up", "Down","Not sig."))
  color_code = list("Up" = "lightgoldenrod",
                    "Down" = "steelblue1",
                    "Not sig." = "gray")
  
  dxr_df$gene_name = gtf_table[match(dxr_df$groupID,gtf_table$gene_id),"gene_name"]
  rownames(dxr_df) = paste(dxr_df$gene_name,dxr_df$featureID,sep = ":")

  volcano_plot_dex <- ggplot(dxr_df,mapping = aes(x = log_2_fold_Change, y=-log10(padj),color = Significance_reg)) +
    geom_point(size=1.75) + 
    xlab("Log 2 fold change") +
    ylab("-log10( p-adjusted )") +
    geom_text(data = dxr_df %>% filter(label == TRUE) , 
              mapping = aes(label = rownames(dxr_df)[which(dxr_df$label == TRUE)]), size = 6, nudge_y = 0.4, check_overlap = T) + 
    
    scale_color_manual(values = color_code) +
    theme(
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      #panel.grid.major = element_blank(), # get rid of major grid
      #panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      #legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
      legend.title = element_text(size = 20, color = "white"),
      legend.key = element_rect(colour = "transparent", fill = "transparent"),
      legend.text = element_text(size = 20, color = "white"),
      axis.line = element_line(color = "white"),
      axis.text = element_text(angle = 45, hjust = 1, size = 17, color = "white"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "white"),
      axis.title = element_text(size = 23, color = "white"))
  output_list$dxr_df = dxr_df
  output_list$volcano_plot = volcano_plot_dex 
  return(output_list)
}




