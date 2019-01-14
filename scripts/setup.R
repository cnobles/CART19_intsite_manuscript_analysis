# setup.R - Setup an appropriate environment for processing the CART19 
# integration site data. This script will check for required packages
# and install those that are missing from a cran mirror (mran).

# Options
options(stringsAsFactors = FALSE, scipen = 99, width = 999)

# Check R-version
Rv <- as.numeric(
  gsub(".", "", paste0(R.Version()$major, R.Version()$minor), fixed = TRUE)
)

if( Rv < 340){
  stop("\n  **Insufficient version of R detected.**\n",
       "  Please use a minimum R version of 3.4.0.")
}

r_libs <- c(
  "igraph", "Matrix", "parallel", "reldist", "reshape2", "RMySQL", "data.table",
  "vegan", "lubridate", "ggrepel", "scales", "grid", "gridExtra", 
  "RColorBrewer", "knitr", "magrittr", "pander", "BiasedUrn", "foreach", 
  "UpSetR", "tidyverse"
)

bioc_libs <- c(
  "BiocGenerics", "Biostrings", "GenomicRanges", "geneRxCluster", "hiAnnotator",
  "sonicLength", "KEGGREST", "GO.db"
)

all_required <- c(r_libs, bioc_libs, "gintools", "spraphal")

is_installed <- suppressMessages(
  sapply(all_required, require, character.only = TRUE)
)

print(
  data.frame(
    "Loaded" = is_installed,
    "R-Package" = names(is_installed)
  ),
  right = FALSE, 
  row.names = FALSE
)

if( !all(is_installed) ){
  
  # Install required R-packages, downloaded from CRAN Mirror snapshot
  get_r_libs <- r_libs[!r_libs %in% row.names(installed.packages())]
  
  if( length(get_r_libs) > 0 ){
    
    install.packages(
      get_r_libs,
      repos = "https://mran.microsoft.com/snapshot/2018-05-01/",
      dependencies = c("Depends", "Imports")
    ) 
    
  }

  # Install BioConductoR-based packages
  suppressMessages(source("https://bioconductor.org/biocLite.R"))
  get_bioc_libs <- bioc_libs[!bioc_libs %in% row.names(installed.packages())]
  
  if( length(get_bioc_libs) > 0 ){
    
    biocLite(
      get_bioc_libs,
      suppressUpdates = TRUE, 
      ask = FALSE,
      siteRepos = "https://mran.microsoft.com/snapshot/2018-05-01/")
    
  }

  # Install developmental based packages from archive
  if( !is_installed["gintools"] ){
    install.packages(
      "packages/gintools_0.1.1.tar.gz", 
      repos = NULL, 
      type = "source"
    )
  }
  
  if( !is_installed["spraphal"] ){
    install.packages(
      "packages/spraphal_0.1.0.tar.gz", 
      repos = NULL, 
      type = "source"
    )
  }

  # Check for installed packages
  all_required <- c(r_libs, bioc_libs, "gintools", "spraphal")
  
  is_installed <- suppressMessages(
    sapply(all_required, require, character.only = TRUE)
  )

  print(
    data.frame(
      "Loaded" = is_installed,
      "R-Package" = names(is_installed)
    ),
    right = FALSE, 
    row.names = FALSE
  )

  if( !all(is_installed) ){
    stop("\n[", paste(Sys.time()), "] Not all required R-packages have been installed. Check dependencies.")
  }else{
    cat("\n[", paste(Sys.time()), "] All required packages installed from extra sources.")
  }
  
}else{
  
  cat("\n[", paste(Sys.time()), "] All required packages installed from sources.")
  
}
