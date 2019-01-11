#    This source code file is a component of the larger INSPIIRED genomic analysis software package.
#    Copyright (C) 2016 Frederic Bushman
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

suppressPackageStartupMessages(library(argparse, quietly=TRUE))

parser <- ArgumentParser(description="Make epigenetic heatmap for sites from database")
parser$add_argument("sample_gtsp", nargs='?', default='sampleName_GTSP.csv')
parser$add_argument("-c", default="./INSPIIRED.yml", help="path to INSPIIRED configuration file.")
parser$add_argument("-o", "--output_dir", type="character", default="epi_heatmap_output",
    help="output folder where genomic heat maps files will be saved")
parser$add_argument("-r", "--ref_genome", type="character", default="hg38", 
    help="reference genome used for all samples (only hg38, hg18, and mm8 are supported at present)")
parser$add_argument("-t", "--cell_types", type="character", 
    help="file with cell types to use: each cell type on separate line(epi_cell_types.R prints all avalable types)")
parser$add_argument("-f", "--file", type="character", default=NULL,
    help="File with sites. CSV format. Included columns: seqnames, start, end, strand, sampleName, refGenome.")

args <- parser$parse_args()

codeDir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly=FALSE), value=T)))

libs <- c("yaml", "survival", "hotROCs", "colorspace", "hiAnnotator", "plyr", "reshape2",
          "DBI", "RMySQL", "BSgenome", "intSiteRetriever", "dplyr", "devtools", "hotROCs")
loaded <- sapply(libs, library, character.only=TRUE, quietly=TRUE)

# Load configuration file
if (!file.exists(args$c)) stop("the configuration file can not be found.")
config <<- yaml.load_file(args$c)

source(file.path(codeDir, "epigeneticFeatures.R"))
source(file.path(codeDir, "epigeneticHeatMapMaker.R"))
source(file.path(codeDir, "utils.R"))

available_epigenetic_ref_genomes <- c ("hg38", "hg18", "mm8")

referenceGenome <- args$ref_genome
if ( ! referenceGenome %in% available_epigenetic_ref_genomes) {
    message("Epigenetic annotation only available for:")
    message(paste(available_epigenetic_ref_genomes, collapse=" "))
    stop(0)
}
heat_map_result_dir <- args$output_dir 
### annotation <- args$annotation_path
annotation <- config$epigeneticDataDirectory
annotation <- file.path(annotation, referenceGenome)

csvfile <- args$sample_gtsp
if( ! file.exists(csvfile) ) stop(csvfile, "not found")
sampleName_GTSP <- read.csv(csvfile)
stopifnot(all(c("sampleName", "GTSP") %in% colnames(sampleName_GTSP)))
message("\nGenerating report from the following sets")
print(sampleName_GTSP)

stopifnot(file.exists(args$file))
sites <- read.csv(args$file)
# Check for required columns
if(!all(c("seqnames", "strand", "position", "sampleName", "refGenome") %in% names(sites))){
  stop("Lacking required columns in input file. See help.")
}
sites <- distinct(sites, seqnames, strand, position, sampleName, refGenome) %>%
  mutate(
    siteID = 1:n(),
    gender = rep("m", n())) %>%
  rename("chr" = seqnames) %>%
  as.data.frame()

histoneorder_for_heatmap <- epigenetic_features()
if ( ! is.null(args$cell_types)) {
    # use only these features
    historder_subset <- read.csv(args$cell_types, header=FALSE, stringsAsFactors=FALSE)$V1
    stopifnot(all(historder_subset %in% histoneorder_for_heatmap)) 
    histoneorder_for_heatmap <- historder_subset
}
message("epigenetic heatmap for the following cell types:")
cat(histoneorder_for_heatmap, sep='\n')

make_epi_heatmap(sampleName_GTSP, referenceGenome, heat_map_result_dir, 
    sites, annotation, histoneorder_for_heatmap)

