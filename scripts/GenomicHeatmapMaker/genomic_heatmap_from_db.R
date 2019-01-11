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

source(file.path(codeDir, "genomicHeatmapMaker.R"))
source(file.path(codeDir, "utils.R"))

libs <- c("argparse", "DBI", "RMySQL", "dplyr", 'yaml', 'RSQLite')
invisible(sapply(libs, library, character.only=TRUE))

parser <- ArgumentParser(description="Make genomic heatmap for sites from database")
parser$add_argument("sample_gtsp", nargs='?', default='sampleName_GTSP.csv')
parser$add_argument("-c", default="./INSPIIRED.yml", help="path to INSPIIRED configuration file.")
parser$add_argument("-o", "--output_dir", type="character", default="heatmap_output",
    help="output folder where genomic heat maps files will be saved")
parser$add_argument("-r", "--ref_genome", type="character", default="hg18", 
    help="reference genome used for all samples")
parser$add_argument("-s", "--sites_group", type="character", default="intsites_miseq.read", 
    help="which group to use for connection")

args <- parser$parse_args()
args

# Load configuration file
if (!file.exists(args$c)) stop("the configuration file can not be found.")
config <<- yaml.load_file(args$c)


referenceGenome <- args$ref_genome
heat_map_result_dir <- args$output_dir 

loaded_ref_genomes <- c ("hg18", "mm9")
if ( ! referenceGenome %in% loaded_ref_genomes) {
    message("Only following genomes are loaded:")
    message(paste(loaded_ref_genomes, collapse=" "))
    message("Install and add new genomes to genomicHeatmapMaker.R")
    message("and add it to loaded_ref_genomes vector in genomic_heatmap_from_db.R")
    stop(0)
}

csvfile <- args$sample_gtsp
if( ! file.exists(csvfile) ) stop(csvfile, " not found")
sampleName_GTSP <- read.csv(csvfile)
stopifnot(all(c("sampleName", "GTSP") %in% colnames(sampleName_GTSP)))
message("\nGenerating report from the following sets")
print(sampleName_GTSP)

# Connect to the database
if (config$dataBase == 'mysql'){
   stopifnot(file.exists("~/.my.cnf"))
   stopifnot(file.info("~/.my.cnf")$mode == as.octmode("600"))
   dbConn <- dbConnect(MySQL(), group=config$mysqlConnectionGroup)
   info <- dbGetInfo(dbConn)
   connection <- src_sql("mysql", dbConn, info = info)
}else if (config$dataBase == 'sqlite') {
   dbConn <- dbConnect(RSQLite::SQLite(), dbname=config$sqliteIntSitesDB)
   info <- dbGetInfo(dbConn)
   connection <- src_sql("sqlite", dbConn, info = info)
   dbConn2 <- dbConnect(RSQLite::SQLite(), dbname=config$sqliteSampleManagement)
   info2 <- dbGetInfo(dbConn2)
   connection2 <- src_sql("sqlite", dbConn2, info = info2)
} else { stop('Can not establish a connection to the database') }

make_heatmap(sampleName_GTSP, referenceGenome, heat_map_result_dir, connection)
