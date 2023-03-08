s
#!/usr/bin/env Rscript

################################################################################
##                                                                            ##
##                                NANOPOREATA                                 ##
##                                                                            ##
################################################################################

# ______________________________________________________________________________
# LIBRARIES ####
options(repos = list(CRAN="http://cran.rstudio.com/"))
source("install.R", local = T)

# Save list of required R packages + version
package_versions = data.table::rbindlist(lapply(sessionInfo()$otherPkgs, function(i) data.frame("Version" = i$Version)), idcol = "pkg")
write.table(package_versions, "NanopoReaTA_Rpackage_versions.txt", sep = "\t", col.names = T, row.names = F, quote = F)

# ______________________________________________________________________________
# SETTINGS ####
options(shiny.maxRequestSize = 30*1024^2)

# ______________________________________________________________________________
# FUNCTIONS ####
## DEA ####
source("server/R_scripts/dea_function.R", local = TRUE)

## DTU and DTE #### 
source("server/R_scripts/dtu_and_dte_function.R", local = TRUE)

## READ LENGTH DISTRIBUTION ####
source("server/R_scripts/read_length_distribution_plots.R", local = TRUE)

## GENE WISE ANALYSIS ####
source("server/R_scripts/gene_wise_analysis_function.R", local  = TRUE)

## INFER EXPERIMENT ####
source("server/R_scripts/infer_experiment_plots.R", local = TRUE)
#____________________________________________________________________________
# FRONTEND
source("ui/ui.R", local = TRUE)
# ______________________________________________________________________________
# BACKEND
source("server/server_docker.R", local = TRUE)

# ______________________________________________________________________________
# LAUNCH APP
APP <- shinyApp(ui, server)
runApp(APP,host = '0.0.0.0',port=8080)

