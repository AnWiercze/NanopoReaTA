inner_var_plot_per_sample <- function(table = data.frame()){
            if (nrow(table) == 0){
                return(ggplot() + theme_void())
            }
            samples = colnames(table)
            table["Iteration"] = as.numeric(rownames(table)) 

            data = c()

            for (i in samples){
              tmp_data = table["Iteration"]
              tmp_data["Difference_of_ratios"] = table[i]
              tmp_data["Sample"] = i
              data = rbind(data,tmp_data)
            }
            data
            if ((nrow(data) / length(samples)) < 20) {
                p = ggline(data,x="Iteration",y="Difference_of_ratios", color="Sample") +
                    theme(
                        # rect = element_rect(fill = "transparent"),
                         panel.background = element_rect(fill = 'transparent', color = "white"), # bg of the panel
                         plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                         legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                        # legend.box.background = element_rect(fill = "transparent"),
                         legend.title = element_text(size = 16, color = "white"),
                         legend.key = element_rect(colour = "transparent", fill = "transparent"),
                         legend.text = element_text(size = 14, color = "white"),
                         axis.text = element_text(angle = 45, hjust = 1, size = 14, color = "white"),
                         plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "white"),
                         axis.title = element_text(size = 23, color = "white"),
                         axis.line = element_line(color = "white"),
                         axis.ticks = element_line(color = "white")) +
                  scale_x_discrete(limits = as.character(seq(1, 20)))
            } else {
                p = ggline(data,x="Iteration",y="Difference_of_ratios", color="Sample") +
                    theme(
                        # rect = element_rect(fill = "transparent"),
                        panel.background = element_rect(fill = 'transparent', color = "white"), # bg of the panel
                        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                        # legend.box.background = element_rect(fill = "transparent"),
                        legend.title = element_text(size = 16, color = "white"),
                        legend.key = element_rect(colour = "transparent", fill = "transparent"),
                        legend.text = element_text(size = 14, color = "white"),
                        axis.text = element_text(angle = 45, hjust = 1, size = 14, color = "white"),
                        plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "white"),
                        axis.title = element_text(size = 23, color = "white"),
                        axis.line = element_line(color = "white"),
                        axis.ticks = element_line(color = "white"))
            }
            p = ggpar(p,legend = "top")
            return(p)
}         

inner_var_plot_per_condition <- function(table = data.frame(), metadata_table, colors=c("#00AFBB", "#E7B800")){
            if (nrow(table) == 0){
                return(ggplot() + theme_void())
            }
            samples = colnames(table)
            table["Iteration"] = as.numeric(rownames(table)) 

            data = c()
            conditions = unique(metadata_table$Condition)

            for (i in conditions){
              tmp_metadata = metadata_table[which(metadata_table$Condition == i),]
              sum = c(0)
              count = 0
              for (j in tmp_metadata$Samples){
                  sum = sum + table[j]
                  count = count + 1
              }
              condition_values = sum/c(count)
              tmp_data = table["Iteration"]
              tmp_data["Difference_of_ratios"] = condition_values 
              tmp_data["Condition"] = i
              data = rbind(data,tmp_data)
            }
            
            if ((nrow(data) / length(tmp_metadata$Samples)) < 20){
                p = ggline(data,x="Iteration",y="Difference_of_ratios", color="Condition", palette=colors) +
                    theme(
                        panel.background = element_rect(fill = "transparent"), # bg of the panel
                        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                        legend.title = element_text(size = 16, color = "white"),
                        legend.key = element_rect(colour = "transparent", fill = "transparent"),
                        legend.text = element_text(size = 14, color = "white"),
                        axis.text = element_text(angle = 45, hjust = 1, size = 14, color = "white"),
                        plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "white"),
                        axis.title = element_text(size = 23, color = "white"),
                        axis.line = element_line(color = "white"),
                        axis.ticks = element_line(color = "white")) +
                    scale_x_discrete(limits = as.character(seq(1, 20)))
            } else {
                p = ggline(data,x="Iteration",y="Difference_of_ratios", color="Condition", palette=colors) +
                    theme(
                        panel.background = element_rect(fill = "transparent"), # bg of the panel
                        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                        legend.title = element_text(size = 16, color = "white"),
                        legend.key = element_rect(colour = "transparent", fill = "transparent"),
                        legend.text = element_text(size = 14, color = "white"),
                        axis.text = element_text(angle = 45, hjust = 1, size = 14, color = "white"),
                        plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "white"),
                        axis.title = element_text(size = 23, color = "white"),
                        axis.line = element_line(color = "white"),
                        axis.ticks = element_line(color = "white")) 
            }
            p = ggpar(p,legend = "top")
            return(p)
}         


total_genes_counted_plot_per_sample <- function(table = data.frame()){
            if (nrow(table) == 0){
                return(ggplot() + theme_void())
            }
            samples = colnames(table)
            table["Iteration"] = as.numeric(rownames(table)) 

            data = c()
            
            for (i in samples){
              tmp_data = table["Iteration"]
              # print(table[i])
              tmp_data["Counts"] = table[i]
              tmp_data["Sample"] = i
              data = rbind(data,tmp_data)
            }
            if ((nrow(data) / length(samples)) < 20){
            q = ggline(data,x="Iteration",y="Counts", color="Sample") + 
                theme(
                    panel.background = element_rect(fill = "transparent"), # bg of the panel
                    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                    legend.title = element_text(size = 16, color = "white"),
                    legend.key = element_rect(colour = "transparent", fill = "transparent"),
                    legend.text = element_text(size = 14, color = "white"),
                    axis.text = element_text(angle = 45, hjust = 1, size = 14, color = "white"),
                    plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "white"),
                    axis.title = element_text(size = 23, color = "white"),
                    axis.line = element_line(color = "white"),
                    axis.ticks = element_line(color = "white")) + 
                scale_x_discrete(limits = as.character(seq(1, 20)))
            } else {
                q = ggline(data,x="Iteration",y="Counts", color="Sample") + 
                    theme(
                        panel.background = element_rect(fill = "transparent"), # bg of the panel
                        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                        legend.title = element_text(size = 16, color = "white"),
                        legend.key = element_rect(colour = "transparent", fill = "transparent"),
                        legend.text = element_text(size = 14, color = "white"),
                        axis.text = element_text(angle = 45, hjust = 1, size = 14, color = "white"),
                        plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "white"),
                        axis.title = element_text(size = 23, color = "white"),
                        axis.line = element_line(color = "white"),
                        axis.ticks = element_line(color = "white"))
            }
            q = ggpar(q,legend = "top")
            return(q)
}


total_genes_counted_plot_per_condition <- function(table = data.frame(), metadata_table, colors=c("#00AFBB", "#E7B800")){
            if (nrow(table) == 0){
                return(ggplot() + theme_void())
            }
            samples = colnames(table)
            table["Iteration"] = as.numeric(rownames(table)) 
            data = c()
            conditions = unique(metadata_table$Condition)

            for (i in conditions){
              tmp_metadata = metadata_table[which(metadata_table$Condition == i),]
              sum = c(0)
              #print(sum)
              count = 0
              #print(tmp_metadata$Samples)
              for (j in tmp_metadata$Samples){
                #print("table J")
                #print(table[j])
                sum = sum + table[j]
                count = count + 1
              }
              # print(sum)
              condition_values = sum/c(count)
              # print(condition_values)
              tmp_data = table["Iteration"]
              tmp_data["Counts"] = condition_values 
              tmp_data["Condition"] = i
              #print(tmp_data)
              data = rbind(data,tmp_data)
            }
            # print(data)
            if ((nrow(data) / length(tmp_metadata$Samples)) < 20){
                p = ggline(data,x="Iteration",y="Counts", color="Condition", palette=colors) + 
                      theme(
                        panel.background = element_rect(fill = "transparent"), # bg of the panel
                        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                        legend.title = element_text(size = 16, color = "white"),
                        legend.key = element_rect(colour = "transparent", fill = "transparent"),
                        legend.text = element_text(size = 14, color = "white"),
                        axis.text = element_text(angle = 45, hjust = 1, size = 14, color = "white"),
                        plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "white"),
                        axis.title = element_text(size = 23, color = "white"),
                        axis.line = element_line(color = "white"),
                        axis.ticks = element_line(color = "white")) +
                    scale_x_discrete(limits = as.character(seq(1, 20)))
            } else {
                p = ggline(data,x="Iteration",y="Counts", color="Condition", palette=colors) + 
                    theme(
                        panel.background = element_rect(fill = "transparent"), # bg of the panel
                        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                        legend.title = element_text(size = 16, color = "white"),
                        legend.key = element_rect(colour = "transparent", fill = "transparent"),
                        legend.text = element_text(size = 14, color = "white"),
                        axis.text = element_text(angle = 45, hjust = 1, size = 14, color = "white"),
                        plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "white"),
                        axis.title = element_text(size = 23, color = "white"),
                        axis.line = element_line(color = "white"),
                        axis.ticks = element_line(color = "white")
                        )
            }
            p = ggpar(p,legend = "top")
            return(p)
}         


inner_var_plot_per_sample.download <- function(table = data.frame()){
    if (nrow(table) == 0){
        return(ggplot() + theme_void())
    }
    samples = colnames(table)
    table["Iteration"] = as.numeric(rownames(table)) 
    
    data = c()
    
    for (i in samples){
        tmp_data = table["Iteration"]
        tmp_data["Difference_of_ratios"] = table[i]
        tmp_data["Sample"] = i
        data = rbind(data,tmp_data)
    }
    data
    p = ggline(data,x="Iteration",y="Difference_of_ratios", color="Sample") +
        theme(
            legend.title = element_text(size = 16, color = "black"),
            legend.text = element_text(size = 14, color = "black"),
            axis.text = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "black"),
            axis.title = element_text(size = 23, color = "black"),
            axis.line = element_line(color = "black"),
            axis.ticks = element_line(color = "black"))
    p = ggpar(p,legend = "top")
    return(p)
}         

inner_var_plot_per_condition.download <- function(table = data.frame(), metadata_table, colors){
    if (nrow(table) == 0){
        return(ggplot() + theme_void())
    }
    samples = colnames(table)
    table["Iteration"] = as.numeric(rownames(table)) 

    data = c()
    conditions = unique(metadata_table$Condition)

    for (i in conditions){
        tmp_metadata = metadata_table[which(metadata_table$Condition == i),]
        sum = c(0)
        count = 0
        for (j in tmp_metadata$Samples){
            sum = sum + table[j]
            count = count + 1
        }
        condition_values = sum/c(count)
        tmp_data = table["Iteration"]
        tmp_data["Difference_of_ratios"] = condition_values 
        tmp_data["Condition"] = i
        data = rbind(data,tmp_data)
    }
    
    p = ggline(data,x="Iteration",y="Difference_of_ratios", color="Condition", palette=colors) +
        theme(
             legend.title = element_text(size = 16, color = "black"),
            legend.text = element_text(size = 14, color = "black"),
            axis.text = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "black"),
            axis.title = element_text(size = 23, color = "black"),
            axis.line = element_line(color = "black"),
            axis.ticks = element_line(color = "black"))
    p = ggpar(p,legend = "top")
    return(p)
}         


total_genes_counted_plot_per_sample.download <- function(table = data.frame()){
    if (nrow(table) == 0){
        return(ggplot() + theme_void())
    }
    samples = colnames(table)
    table["Iteration"] = as.numeric(rownames(table)) 
    
    data = c()
    
    for (i in samples){
        tmp_data = table["Iteration"]
        tmp_data["Counts"] = table[i]
        tmp_data["Sample"] = i
        data = rbind(data,tmp_data)
    }
    
    q = ggline(data,x="Iteration",y="Counts", color="Sample") + 
        theme(
            legend.title = element_text(size = 16, color = "black"),
            legend.text = element_text(size = 14, color = "black"),
            axis.text = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "black"),
            axis.title = element_text(size = 23, color = "black"),
            axis.line = element_line(color = "black"),
            axis.ticks = element_line(color = "black"))
    q = ggpar(q,legend = "top")
    return(q)
}


total_genes_counted_plot_per_condition.download <- function(table = data.frame(), metadata_table, colors){
    if (nrow(table) == 0){
        return(ggplot() + theme_void())
    }
    samples = colnames(table)
    table["Iteration"] = as.numeric(rownames(table)) 
    data = c()
    conditions = unique(metadata_table$Condition)

    for (i in conditions){
        tmp_metadata = metadata_table[which(metadata_table$Condition == i),]
        sum = c(0)
        count = 0
        for (j in tmp_metadata$Samples){
            sum = sum + table[j]
            count = count + 1
        }
        condition_values = sum/c(count)
        tmp_data = table["Iteration"]
        tmp_data["Counts"] = condition_values 
        tmp_data["Condition"] = i
        data = rbind(data,tmp_data)
    }
    
    p = ggline(data,x="Iteration",y="Counts", color="Condition", palette=colors) + 
        theme(
            legend.title = element_text(size = 16, color = "black"),
            legend.text = element_text(size = 14, color = "black"),
            axis.text = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "black"),
            axis.title = element_text(size = 23, color = "black"),
            axis.line = element_line(color = "black"),
            axis.ticks = element_line(color = "black"))
    p = ggpar(p,legend = "top")
    return(p)
}         

