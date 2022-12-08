# source of this script:
#    Title: NanopoReaTA
#    Author: all
#    Date: 16.08.2022
#    Availability: https://github.com/AnWiercze/testApp


################################################################################
# Check that the currently-installed version of R
# is the correct version
################################################################################

options(repos = list(CRAN="http://cran.rstudio.com/"))

R_min_version = "4.1.2"
R_version = paste0(R.Version()$major, ".", R.Version()$minor)
if(compareVersion(R_version, R_min_version) == -1){
 stop("You need to have at least version 4.1.3 of R to run the app.\n", 
      "Launch should fail.\n",
      "Go to http://cran.r-project.org/ and install version 4.1.3 of R or higher.")
}

#Check if BiocManager is installed and install otherwise
availpacks = .packages(all.available = TRUE)
if (!("BiocManager" %in% availpacks)){
  install.packages("BiocManager")
  
}

################################################################################
# Install basic required packages if not available/installed.
################################################################################
install_missing_packages = function(pkg, version = NULL, verbose = TRUE){
  availpacks = .packages(all.available = TRUE)
  require("BiocManager")
  missingPackage = FALSE
  if(!any(pkg %in% availpacks)){
    if(verbose){
      message("The following package is missing.\n",
              pkg, "\n",
              "Installation will be attempted...")
    }
    missingPackage <- TRUE
  }
  if(!is.null(version) & !missingPackage){
    # version provided and package not missing, so compare.
    if( compareVersion(a = as.character(packageVersion(pkg)),
                       b = version) < 0 ){
      if(verbose){
        message("Current version of package\n", 
                pkg, "\t", 
                packageVersion(pkg), "\n",
                "is less than required.
                Update will be attempted.")
      }
      missingPackage <- TRUE
    }
  }
  if(missingPackage){
      BiocManager::install(pkg, update=FALSE)
  }
}
################################################################################
# Define list of package names and required versions.
################################################################################
deppkgs = readRDS("NanopoReaTA_Rpackages.RDS")
#deppkgs = readRDS("tryout_versions.RDS")

# Loop on package check, install, update
pkg1 = mapply(install_missing_packages,
              pkg = names(deppkgs), 
              version = deppkgs,
              MoreArgs = list(verbose = TRUE), 
              SIMPLIFY = FALSE,
              USE.NAMES = TRUE)
################################################################################
# Load packages 
################################################################################
for(i in names(deppkgs)){
  suppressPackageStartupMessages({library(i, character.only = TRUE)})
  message(i, " package version:\n", packageVersion(i))
}
################################################################################
