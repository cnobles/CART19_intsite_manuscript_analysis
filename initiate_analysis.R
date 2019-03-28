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

# In our experience, we've identified occational sites that are known for
# mispriming, their position identifiers are listed here and are excluded from
# the analysis
potential_mispriming <- c(
  "chr1-111167195", "chr1-17174312", "chr1-199385510", "chr1-20336952", 
  "chr1-222611455", "chr1-232153262", "chr1-246907579", "chr1-32161842", 
  "chr1-62345880", "chr1-88021817", "chr1-98271696", "chr1+163797078", 
  "chr1+170406779", "chr1+196997028", "chr1+230477392", "chr1+37108459", 
  "chr1+46670094", "chr1+54742232", "chr1+65480833", "chr1+91362383", 
  "chr1+92317581", "chr2-138790717", "chr2-187281821", "chr2-213462988", 
  "chr2-22258798", "chr2-224307606", "chr2-225008345", "chr2-80960242", 
  "chr2-98492415", "chr2+124136527", "chr2+137430498", "chr2+14869402", 
  "chr2+152041816", "chr2+16347997", "chr2+180744900", "chr2+188290705", 
  "chr2+232108747", "chr2+433020", "chr2+7274721", "chr2+80735962", 
  "chr2+9650737", "chr3-176801222", "chr3-32554245", "chr3-47305137", 
  "chr3-56858480", "chr3-66663178", "chr3+128337787", "chr3+14036591", 
  "chr3+171396716", "chr3+171497564", "chr3+17869203", "chr3+184884914", 
  "chr3+193697301", "chr3+21594102", "chr3+3244449", "chr3+73437766", 
  "chr3+76579041", "chr3+80583046", "chr3+97455063", "chr4-109191857", 
  "chr4-156033113", "chr4-171835401", "chr4-183829", "chr4-20562862", 
  "chr4-23756717", "chr4-35878780", "chr4-54890595", "chr4-595582", 
  "chr4-66697715", "chr4+141075796", "chr4+37791681", "chr4+43501444", 
  "chr5-113738443", "chr5-151311821", "chr5-29684662", "chr5-51240157", 
  "chr5-71561837", "chr5-81202038", "chr5+122942607", "chr5+126064143", 
  "chr5+133355340", "chr5+144653533", "chr5+26071689", "chr5+41977260", 
  "chr5+60346535", "chr5+81591158", "chr6-112403203", "chr6-17678597", 
  "chr6-19603897", "chr6-51795350", "chr6+114989283", "chr6+135784714", 
  "chr6+136442896", "chr6+18418179", "chr6+24182787", "chr6+34364358", 
  "chr6+47058858", "chr6+49058179", "chr6+79138366", "chr7-103762013", 
  "chr7-106141866", "chr7-116985749", "chr7-126568788", "chr7-139261550", 
  "chr7-142174493", "chr7-151890389", "chr7-16528708", "chr7-4209700", 
  "chr7-71655708", "chr7-7701466", "chr7-77856029", "chr7+10384476", 
  "chr7+125751530", "chr7+143244258", "chr7+14992752", "chr7+150525009", 
  "chr7+21002952", "chr7+67681333", "chr7+87013270", "chr7+96131624", 
  "chr8-124129504", "chr8-131186001", "chr8-37102591", "chr8-82308317", 
  "chr8-90063638", "chr8-95000415", "chr8+100540408", "chr8+108712491", 
  "chr8+141148698", "chr8+3311190", "chr8+60387404", "chr8+65074905", 
  "chr9-109489420", "chr9-36018607", "chr9-7602396", "chr9-82101389", 
  "chr9+102255244", "chr9+118088336", "chr9+75214874", "chr9+86865891", 
  "chr9+95003739", "chr10-103386926", "chr10-115649405", "chr10-27180287", 
  "chr10-5035518", "chr10-84444163", "chr10+58362781", "chr10+77681314", 
  "chr10+88967386", "chr10+935660", "chr10+99029060", "chr11-120463032", 
  "chr11-123354241", "chr11-30268387", "chr11-59789938", "chr11-69462663", 
  "chr11+18624319", "chr11+24773832", "chr11+27750570", "chr11+29110926", 
  "chr11+36104092", "chr11+87266909", "chr11+88642457", "chr11+96431373", 
  "chr12-112101436", "chr12-39548196", "chr12-41002459", "chr12-91551285", 
  "chr12+29335943", "chr12+4280324", "chr12+4876938", "chr12+98603313", 
  "chr13-24256934", "chr13-38737440", "chr13-52853795", "chr13-58786705", 
  "chr13+22014042", "chr13+25972050", "chr13+26395024", "chr13+95919557", 
  "chr14-38012095", "chr14-85146958", "chr14+104956207", "chr14+50214279", 
  "chr14+84816180", "chr14+98324668", "chr15-31567573", "chr15-42179737", 
  "chr15-68005864", "chr15+71312727", "chr16-35905836", "chr16-81876928", 
  "chr16+2230954", "chr17-31227950", "chr17-44041761", "chr17-73686205", 
  "chr17-80568713", "chr17+14319231", "chr17+37150381", "chr17+40640172", 
  "chr17+65038839", "chr17+7656638", "chr18-32708283", "chr18-36634891", 
  "chr18-40776203", "chr18-8829148", "chr18+27211376", "chr18+31008461", 
  "chr18+4450406", "chr19-30847912", "chr19-40334020", "chr19-57808512", 
  "chr19+14381221", "chr19+1492075", "chr19+23671602", "chr20-16215014", 
  "chr20-3753172", "chr20-4223765", "chr20-50662792", "chr20-63204711", 
  "chr20+16798194", "chr20+31302607", "chr20+58550365", "chr20+62951002", 
  "chr21+29891265", "chr21+42691561", "chr22-33610670", "chr22-37189964", 
  "chr22-37323913", "chr22-42819667", "chrX-117238952", "chrX-127607004", 
  "chrX-145868827", "chrX-21400975", "chrX-22629243", "chrX-30400966", 
  "chrX-32704035", "chrX-52085042", "chrX-6788253", "chrX+145623876", 
  "chrX+36593140", "chrX+67072314", "chrX+71385589", "chrX+88391839", 
  "chrY+15068888"
)

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
      "-r hg38 -c utils/INSPIIRED.yml ",
      "-u utils"
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