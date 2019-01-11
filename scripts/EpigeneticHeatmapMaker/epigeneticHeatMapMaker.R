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

make_epi_heatmap <- function(sampleName_GTSP, referenceGenome, output_dir, 
    connection, annotPath, histoneorder
) {
    if(class(connection) == "character"){
      sites <- get_sites_controls_from_db(sampleName_GTSP, referenceGenome, connection)
    }else{
      sites <- get_sites_controls_from_file(sampleName_GTSP, referenceGenome, connection)
    }
    sites_to_epigenetic_heatmap(sites, referenceGenome, sampleName_GTSP, output_dir, annotPath, histoneorder)
}

sites_to_epigenetic_heatmap <- function(sites, referenceGenome, sampleName_GTSP, output_dir, annotPath, histoneorder) {
    #### set parameters for annotations & analysis ####
    if( ! file.exists(annotPath)) {
        stop("Required annotation directory doesn't exist: ", annotPath)
    }
    inoutNuc <- 'yes'
    windows <- c(10000)
    #### See if combinations of setName-epigeneticFactor-window exist ####
    ## get counts of all sites in samples ##
    todocombos <- expand.grid(setName=unique(sites$sampleName), histones=histoneorder, 
                              windows=windows, stringsAsFactors=FALSE)
    todocombos$histone_window <- with(todocombos, 
                                      paste(histones, getWindowLabel(windows),sep="."))
    todocombos$todo <- TRUE
    rows <- todocombos$todo
    todoHistones <- unique(todocombos$histones[rows])
    todoWindows <- unique(todocombos$windows[rows])
    if(length(sites)>0) {
        for(f in todoHistones) {
            message(f)
            load(file.path(annotPath, paste0(f, ".RData")))
            if(length(sites)>1e6) {
                sites <- getFeatureCounts(sites, epigenData, f, width=todoWindows,
                                          doInChunks=TRUE)
            } else {
                sites <- getFeatureCounts(sites, epigenData, f, width=todoWindows)
            }
            if(tolower(inoutNuc)=="yes" & referenceGenome=="hg18" & 
                   f %in% c("ActivatedNucleosomes","RestingNucleosomes")) {        
                methinout <- paste(f,"inout",sep=".")
                epigenData$name <- "bore"
                sites <- getSitesInFeature(sites, epigenData, methinout, asBool=TRUE)
            }
            rm(epigenData)
        }
    }

    if (config$rocControls == 'unmatched')
    {
      sites_to_ROC_ordinary(sites, sampleName_GTSP, output_dir)
    } else if (config$rocControls == 'matched') {
      sites_to_ROC_matched(sites, sampleName_GTSP, output_dir)
    } else {
      stop('Error, rocControls is not properly defined in the configuration file.')
    }
}
