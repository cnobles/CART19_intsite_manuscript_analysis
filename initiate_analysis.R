# initiate_analysis.R
# usage: Rscript initiate_analysis.R -c <cores>
# params: -c NUM    Number of CPU cores to used during processing, default 30.
#                   Will not use more cores than are available on the system.
#                   
# Description: To initiate the CART19 analysis, please execute this script 
# (using the Rscript command) in the repository directory. 
# (i.e. ~/CART19_intsite_manuscript_analysis)

# Setup ----
cat("[", paste(Sys.time()), "] Starting setup and checking installed packages.\n")

## Working paths:
workingDir <- getwd()
scriptDir <- file.path(workingDir, "scripts")
utilsDir <- file.path(workingDir, "utils")
outputDir <- file.path(workingDir, "data")
vcnDir <- file.path(workingDir, "data", "vcn_data")

which_Rscript <- "/home/opt/R-3.4.0/bin/Rscript"
if( !file.exists(which_Rscript) ) which_Rscript <- "Rscript"

if( !dir.exists(outputDir) ){
  system(paste0("mkdir ", outputDir))
}

source(file.path(scriptDir, "setup.R"))

# Initialize processing from unprocessed integration sites
numCores <- as.numeric(commandArgs(trailingOnly = TRUE)[
  grep("-c", commandArgs(trailingOnly = TRUE)) + 1 
])

if( length(numCores) == 0 ) numCores <- parallel::detectCores()

numCores <- min(numCores, detectCores())

cat("[", paste(Sys.time()), "] Number of cores to use: ", numCores, "\n")

# Starting initial processing
cat("[", paste(Sys.time()), "] Initial processing started...\n")
source(
  file = file.path(scriptDir, "initial_cart19_data_processing.R"), 
  local = new.env()
)

# Initial processing only needs to be run to generate all the data needed for
# analysis. This takes quite a while to run (several hours) due to the amount 
# of data and standardizing algorithms that must be done for quality control.

# As part of the initial processing, contamination filtering only needs to be 
# run once to filter the unique sites. Once the output file is saved, it does 
# not need to be reprocessed.

# The following scripts will process the data into a variety of forms, and save
# the output to the data directory.

# Create a processing environment

# Setup processing environment ----
cat("[", paste(Sys.time()), "] Running 'proc_env_setup.R'... (Processing step 1 of 7)\n")
suppressMessages(source(file = file.path(scriptDir, "proc_env_setup.R")))

# Condense and annotate unique sites ----
cat("[", paste(Sys.time()), "] Running 'uniq_process_and_annotate.R'... (Processing step 2 of 7)\n")
suppressMessages(source(file = file.path(scriptDir, "uniq_process_and_annotate.R")))

# Process scan statistics on unique sites ----
cat("[", paste(Sys.time()), "] Running 'process_scan_statistics.R'... (Processing step 3 of 7)\n")
suppressMessages(source(file = file.path(scriptDir, "process_scan_statistics.R")))

# Summarize unique integration sites across patient, timepoint, ... ----
cat("[", paste(Sys.time()), "] Running 'summarize_intsites.R'... (Processing step 4 of 7)\n")
suppressMessages(source(file = file.path(scriptDir, "summarise_intsites.R")))

# Identify and consolidate genes of interest from cluster data ----
cat("[", paste(Sys.time()), "] Running 'consolidate_genes_of_interest.R'... (Processing step 5 of 7)\n")
suppressMessages(source(file = file.path(scriptDir, "consolidate_genes_of_interest.R")))

# Determine total genic impact (only sites within genes) ----
cat("[", paste(Sys.time()), "] Running 'genic_impact.R'... (Processing step 6 of 7)\n")
suppressMessages(source(file = file.path(scriptDir, "genic_impact.R")))

# Determine the intergenic impact (only non-gene sites) ----
cat("[", paste(Sys.time()), "] Running 'intergenic_impact.R'... (Processing step 7 of 7)\n")
suppressMessages(source(file = file.path(scriptDir, "intergenic_impact.R")))

# HeatmapMakers ----
# After analytical processing, some data needs to be processed by utilities from
# INSPIIRED. The flowing sets up the files relevent to the heatmap analysis 
# focusing on genomic or epigenetic features.

# Generate heatmaps ----
# Select TDN and Patient sites (d28) for heatmap
# additionally split by responders and non-responders

cat(
  "[", paste(Sys.time()), "] Generating heatmaps using GenomicHeatmapMaker and EpigeneticHeatmapMaker...\n"
)

req_data <- c("tdn_sites", "timepoint_sites")

if( !all(sapply(req_data, exists)) ){
  tdn_sites <- readRDS(file.path(outputDir, "cart19_tdn_sites.rds"))
  timepoint_sites <- readRDS(file.path(outputDir, "cart19_timepoint_sites.rds"))
}
 
std_data <- dplyr::bind_rows(list(
    "TDN" = as.data.frame(tdn_sites, row.names = NULL), 
    "D28" = as.data.frame(timepoint_sites, row.names = NULL)),
    .id = "type"
  ) %>%
  dplyr::select(type, patient, celltype, timepoint, specimen, posid) %>%
  dplyr::filter(type == "D28" & timepoint == "d28" | type == "TDN") %>%
  dplyr::left_join(patient_data, by = "patient") %>%
  dplyr::mutate(
    seqnames = stringr::str_extract(posid, "[\\w]+"),
    strand = stringr::str_extract(posid, "[+-]"),
    pos = as.numeric(stringr::str_extract(posid, "[0-9]+$")),
    setname = paste0(type, " - ", determinant_response)) %>%
  dplyr::select(seqnames, strand, pos, setname, posid) %>% 
  dplyr::filter(seqnames %in% c(paste0("chr", 1:22), "chrX", "chrY")) %>%
  as.data.frame()

std_data <- split(std_data, std_data$setname)
std_heatmap_data <- dplyr::bind_rows(lapply(std_data, function(df){
    if(nrow(df) < 5000){
      return(df)
    }else{
      idx <- sample(
        seq_len(nrow(df)), 5000, replace = FALSE)
      return(df[idx,])
    }
  })) %>%
  dplyr::mutate(sampleName = setname, position = pos, refGenome = "hg38") %>%
  dplyr::select(seqnames, strand, position, sampleName, refGenome)

std_heatmap_sample_info <- std_heatmap_data %>%
  dplyr::group_by(sampleName) %>%
  dplyr::summarise(
    GTSP = unique(sampleName), 
    patient = "Mock") %>%
  dplyr::ungroup() 

write.csv(
  std_heatmap_data, 
  file.path(outputDir, "manuscript_std_heatmap_sites.csv"), 
  row.names = FALSE, 
  quote = FALSE
)

write.csv(
  std_heatmap_sample_info, 
  file.path(outputDir, "manuscript_std_heatmap_samples.csv"), 
  row.names = FALSE, 
  quote = FALSE
)

rm(std_data, std_heatmap_data, std_heatmap_sample_info)

# Call for Genomic Heatmap
if( !file.exists("data/manuscript_genomic_std_heatmap/roc.res.rds") ){
  
  system(
    paste0(
      which_Rscript, 
      " scripts/GenomicHeatmapMaker/genomic_heatmap_from_file.R ",
      "data/manuscript_std_heatmap_samples.csv ",
      "-f data/manuscript_std_heatmap_sites.csv ",
      "-o data/manuscript_genomic_std_heatmap ",
      "-r hg38 -c utils/INSPIIRED.yml"
    ), 
    wait = TRUE
  )
  
  while( !file.exists("data/manuscript_genomic_std_heatmap/roc.res.rds") ){
    Sys.sleep(10)
  }
  
}

# Call for Epigenetic Heatmap 
if( !file.exists("data/manuscript_epi_std_heatmap/roc.res.rds") ){
  
  system(
    paste0(
      which_Rscript,
      " scripts/EpigeneticHeatmapMaker/epi_heatmap_from_file.R ",
      "data/manuscript_std_heatmap_samples.csv ",
      "-f data/manuscript_std_heatmap_sites.csv ",
      "-t scripts/EpigeneticHeatmapMaker/CD4_epi_types.txt ",
      "-o data/manuscript_epi_std_heatmap ",
      "-c utils/INSPIIRED.yml"
    ),
    wait = TRUE
  )
  
  while( !file.exists("data/manuscript_epi_std_heatmap/roc.res.rds") ){
    Sys.sleep(10)
  }
  
}
  
# Report generation ----
# After heatmap generation has completed, the remaining analysis is focused in
# specific reports, either the manuscript analysis which generates the base 
# figures for the paper, or the genes of interest analysis, which generates the
# supplementary report in the manuscript.

if(
  file.exists("reports/cart19_goi_report.archive.v1.pdf") &
  file.exists("reports/cart19_intsite_analysis.archive.v1.pdf")
  ){
  
  cat("[", paste(Sys.time()), "] Reports found in reports directory. Analysis complete.\n")
  
}else{

  cat("[", paste(Sys.time()), "] Starting report generation, first goi report...\n")
  
  # Generate genes of interest report
  goi_report_template_path <- normalizePath(
    file.path(scriptDir, "cart19_goi_report.archive.v1.Rmd")
  )
  
  rmarkdown::render(
    input = goi_report_template_path,
    output_format = "pdf_document", 
    output_file = "cart19_goi_report.archive.v1.pdf",
    output_dir = "reports"
  )
  
  while( !file.exists("reports/cart19_goi_report.archive.v1.pdf") ){
    Sys.sleep(10)
  }
  
  cat("\n[", paste(Sys.time()), "] Second report for manuscript figures...\n")
  
  # Generate manuscript report
  manuscript_report_template_path <- normalizePath(
    file.path(scriptDir, "cart19_intsite_analysis.archive.v1.Rmd")
  )
  
  rmarkdown::render(
    input = manuscript_report_template_path,
    output_format = "pdf_document", 
    output_file = "cart19_intsite_analysis.archive.v1.pdf",
    output_dir = "reports"
  )
  
  
  Sys.sleep(10)
  
  if( 
    file.exists("reports/cart19_goi_report.archive.v1.pdf") &
    file.exists("reports/cart19_intsite_analysis.archive.v1.pdf") 
  ){
    
    cat("[", paste(Sys.time()), "] Reports generated. Analysis has completed.\n")
    q()
    
  }else{
    
    stop("[", paste(Sys.time()), "] Reports not found in 'reports' directory, check for error.\n")
  
  }

}