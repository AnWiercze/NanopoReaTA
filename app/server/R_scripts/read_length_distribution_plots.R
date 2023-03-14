################################################################################
##                          Length distribution                               ##
################################################################################

# This script creates plots of sample- and group-wise length distributions. 

theme_set(theme_light())
theme_update(
  panel.background = element_rect(fill = "transparent"), # bg of the panel
  plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
  legend.background = element_rect(fill = "transparent"), # get rid of legend bg
  legend.title = element_text(size = 20, color = "white"),
  legend.key = element_rect(colour = "transparent", fill = "transparent"),
  legend.text = element_text(size = 20, color = "white"),
  axis.text = element_text(angle = 45, hjust = 1, size = 17, color = "white"),
  plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "white"),
  axis.title = element_text(size = 23, color = "white"),
  axis.line = element_line(color = "white"),
  axis.ticks = element_line(color = "white"))

createLengthPlots <- function(readLengths_df_filt, metadata, conditionCol, conditions, color_conditions){
  # Extract conditions of interest 
  metadata = metadata[metadata[[conditionCol]] %in% conditions, ]
  
  # Join the metadata information with read length file
  readLengths_df_filt = readLengths_df_filt %>% 
    left_join(metadata, by = c("Sample" = "Samples"))
  # Plot sample-wise distribution of all reads
  sampleWise_All = ggplot(readLengths_df_filt, aes(x=Length, color=Sample))  +
    geom_density() +
    scale_x_continuous(labels=scales::comma) +
    xlab("Read length") +
    ylab("Density") +
    ggtitle("Sample-wise length distribution\n(filtered longest 1 % of reads)") +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Plot group-wise distribution of all reads
  groupWise_All = ggplot(readLengths_df_filt, aes(x = Length, color = Condition)) +
    scale_color_manual(values=color_conditions) + 
    geom_density() +
    scale_x_continuous(labels=scales::comma) +
    xlab("Read length") +
    ylab("Density") +
    ggtitle("Condition-wise length distribution\n(filtered longest 1 % of reads)") +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  return(list(sampleWise_All = sampleWise_All, groupWise_All = groupWise_All))
}
samplewise_read_length.download <- function(readLengths_df_filt, metadata, conditionCol, conditions){
  theme_update(legend.title = element_text(size = 20, color = "black"),
               # legend.key = element_rect(colour = "transparent", fill = "transparent"),
               legend.text = element_text(size = 20, color = "black"),
               axis.text = element_text(angle = 45, hjust = 1, size = 17, color = "black"),
               plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "black"),
               axis.title = element_text(size = 23, color = "black"),
               axis.line = element_line(color = "black"),
               axis.ticks = element_line(color = "black"), 
               panel.background = element_rect(fill = "white"), # bg of the panel
               plot.background = element_rect(fill = "white"), # bg of the plot
               legend.background = element_rect(fill = "white"))
  metadata = metadata[metadata[[conditionCol]] %in% conditions, ]
  
  # Join the metadata information with read length file
  readLengths_df_filt = readLengths_df_filt %>% 
    left_join(metadata, by = c("Sample" = "Samples"))
  
  # Plot sample-wise distribution of all reads
  sampleWise_All = ggplot(readLengths_df_filt, aes(x = Length, color = Sample)) +
    geom_density() +
    scale_x_continuous(labels=scales::comma) +
    xlab("Read length") +
    ylab("Density") +
    ggtitle("All reads: Sample-wise length distribution\n(filtered longest 1 % of reads)") +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(sampleWise_All)
  
}
groupwise_read_length.download <- function(readLengths_df_filt, metadata, conditionCol, conditions, color_conditions){
  metadata = metadata[metadata[[conditionCol]] %in% conditions, ]
  theme_update(legend.title = element_text(size = 20, color = "black"),
               # legend.key = element_rect(colour = "transparent", fill = "transparent"),
               legend.text = element_text(size = 20, color = "black"),
               axis.text = element_text(angle = 45, hjust = 1, size = 17, color = "black"),
               plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "black"),
               axis.title = element_text(size = 23, color = "black"),
               axis.line = element_line(color = "black"),
               axis.ticks = element_line(color = "black"), 
               panel.background = element_rect(fill = "white"), # bg of the panel
               plot.background = element_rect(fill = "white"), # bg of the plot
               legend.background = element_rect(fill = "white"))
  # Join the metadata information with read length file
  readLengths_df_filt = readLengths_df_filt %>% 
    left_join(metadata, by = c("Sample" = "Samples"))
  
  # Plot group-wise distribution of all reads
  groupWise_All = ggplot(readLengths_df_filt, aes(x = Length, color = Condition)) +
    scale_color_manual(values=color_conditions) + 
    geom_density() +
    scale_x_continuous(labels=scales::comma) +
    xlab("Read length") +
    ylab("Density") +
    ggtitle("All reads: Condition-wise length distribution\n(filtered longest 1 % of reads)") +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(groupWise_All)
}


