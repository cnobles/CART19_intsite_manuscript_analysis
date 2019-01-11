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
#

codeDir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly=FALSE), value=T)))

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(description="Make genomic heatmap for sites from database")
parser$add_argument("sample_gtsp", nargs='?', default='sampleName_GTSP.csv')
parser$add_argument("-c", default="./INSPIIRED.yml", help="path to INSPIIRED configuration file.")
parser$add_argument("-o", "--output_dir", type="character", default="heatmap_output",
    help="output folder where genomic heat maps files will be saved")
parser$add_argument("-r", "--ref_genome", type="character", default="hg18", 
    help="reference genome used for all samples")
parser$add_argument("-f", "--file", type="character", default=NULL,
    help="File with sites. CSV format. Included columns: seqnames, strand, position, sampleName, refGenome.")

args <- parser$parse_args()
args

suppressMessages(source(file.path(codeDir, "genomicHeatmapMaker.R")))
suppressMessages(source(file.path(codeDir, "utils.R")))

libs <- c("DBI", "RMySQL", "dplyr", 'yaml', 'RSQLite', 'hotROCs')
null <- suppressMessages(sapply(libs, library, character.only=TRUE))

# Load configuration file
if (!file.exists(args$c)) stop("the configuration file can not be found.")
config <<- yaml.load_file(args$c)

heat_map_result_dir <- args$output_dir 

genome <- grep(args$ref_genome, BSgenome::installed.genomes(), value = TRUE)
  if(length(genome) == 0){
    print("Installed genomes include:")
    print(BSgenome::installed.genomes())
    stop("Selected reference genome not in list.")
  }else if(length(genome) > 1){
    print("Installed genomes include:")
    print(BSgenome::installed.genomes())
    stop(
      "Please be more specific about reference genome. Multiple matches to input.")
  }
suppressMessages(library(genome, character.only = TRUE))
referenceGenome <- args$ref_genome

csvfile <- args$sample_gtsp
if( ! file.exists(csvfile) ) stop(csvfile, " not found.")
sampleName_GTSP <- read.csv(csvfile)
stopifnot(all(c("sampleName", "GTSP") %in% colnames(sampleName_GTSP)))
message("\nGenerating report from the following sets")
print(sampleName_GTSP)

stopifnot(file.exists(args$file))
sites <- read.csv(args$file)
# Check for required columns
if(!all(
  c("seqnames", "strand", "position", "sampleName", "refGenome") %in% 
  names(sites))){
    stop("Lacking required columns in input file. See help.")
}
sites <- distinct(sites, seqnames,strand, position, sampleName, refGenome) %>%
  mutate(
    siteID = 1:n(),
    gender = rep("m", n())) %>%
  rename("chr" = seqnames) %>%
  as.data.frame()

make_heatmap(sampleName_GTSP, referenceGenome, heat_map_result_dir, sites)
