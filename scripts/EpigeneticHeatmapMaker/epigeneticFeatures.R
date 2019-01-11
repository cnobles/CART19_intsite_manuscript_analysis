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

epigenetic_features <- function() {
    # datasets
    meth.sets <- c('H2BK5me1','H3R2me1','H3R2me2','H3K4me1','H3K4me2','H3K4me3','H3K9me1','H3K9me2','H3K9me3','H3K27me1','H3K27me2','H3K27me3','H3K36me1','H3K36me3','H3K79me1','H3K79me2','H3K79me3','H4R3me2','H4K20me1','H4K20me3','PolII','H2AZ','CTCF','H3K4me1-HeLaS3_Untr','H3K4me1-HeLaS3_INFS','H3K4me3-HeLaS3_Untr','H3K4me3-HeLaS3_INFS','H3K27me3-CD133','H3K4me1-CD133','H3K9me1-CD133','H4K20me1-CD133','H3K27me1-CD133','H3K36me3-CD133','H3K4me3-CD133','H3K9me3-CD133','PolII-CD133','H2AZ-CD133','H3K9me2-HeLa','H3K9me3-HeLa','H3K27me3-HeLa','H3K36me3-HeLa','H4K20me3-HeLa','CTCF-HeLa','H4K20me1-HeLa','H2BK5me1-HeLa','H3K4me1-HeLa','H2AZ_hsalt-HeLa','H2AZ_lsalt-HeLa','H3K4me3-HeLa')
    act.sets <- c('H2AK5ac','H2AK9ac','H2BK120ac','H2BK12ac','H2BK20ac','H2BK5ac','H3K14ac','H3K18ac','H3K23ac','H3K27ac','H3K36ac','H3K4ac','H3K9ac','H4K12ac','H4K16ac','H4K5ac','H4K8ac','H4K91ac','H3K27ac-HeLa') 
    tf.sets <- c('NRSF','STAT1-HeLaS3_Untr','STAT1-HeLaS3_INFS','STAT1-HeLaS3_INFS-2','ActivatedNucleosomes','RestingNucleosomes','Rest-CD4-HDAC1','Rest-CD4-HDAC2','Rest-CD4-HDAC3','Rest-CD4-HDAC6','Act-CD4-HDAC6','Rest-CD4-CBP','Rest-CD4-PCAF','Rest-CD4-p300','Rest-CD4-MOF','Rest-CD4-Tip60','Act-CD4-Tip60','H3_3-HeLa','p300-HeLa')

    #species <- switch(referenceGenome, hg18="Homo sapiens", mm8="Mus musculus")
    #sites <- keepStandardChromosomes(sites, style="UCSC", species=species)

    #### Set the order of histone modifications to be displayed on the heatmap ####
    display.order <- c("H3K9me2", "H3K9me2-HeLa", "H3K9me3", "H3K9me3-MEF", 
                       "H3K9me3-CD133", "H3K9me3-HeLa", "H3K27me2", "H3K27me3", 
                       "H3K27me3-MEF", "H3K27me3-CD133", "H3K27me3-HeLa", "H3K14ac",
                       "H2AK9ac", "H3K23ac", "H3K36me1", "H3K36me3", "H3K36me3-MEF",
                       "H3K36me3-CD133", "H3K36me3-HeLa", "H3K27me1", "H3K27me1-CD133",
                       "H3R2me1", "H3R2me2", "H4R3me2", "H4K20me3", "H4K20me3-HeLa",
                       "H2AK5ac", "H4K16ac", "H4K12ac", "CTCF", "CTCF-HeLa", 
                       "CTCF-Jurkat", "PolII", "PolII-CD133", "H3K79me1", "H3K79me2",
                       "H3K79me3", "H4K20me1", "H4K20me1-CD133", "H4K20me1-HeLa", 
                       "H2BK5me1", "H2BK5me1-HeLa", "H3K4me1", "H3K4me1-CD133", 
                       "H3K4me1-HeLa", "H3K4me1-HeLaS3_Untr", "H3K4me1-HeLaS3_INFS",
                       "H3K4me2", "H3K9me1", "H3K9me1-CD133", "H4K8ac", "H4K5ac", 
                       "H2AZ", "H2AZ-CD133", "H2AZ_hsalt-HeLa", "H2AZ_lsalt-HeLa", 
                       "H3K4me3", "H3K4me3-MEF", "H3K4me3-CD133", "H3K4me3-HeLa",
                       "H3K4me3-HeLaS3_Untr", "H3K4me3-HeLaS3_INFS", "H2BK12ac", 
                       "H3K9ac", "H3K18ac", "H3K27ac", "H3K27ac-HeLa", "H2BK5ac", 
                       "H3K36ac", "H3K4ac", "H2BK20ac", "H4K91ac", "H2BK120ac", 
                       "NRSF", "STAT1-HeLaS3_Untr", "STAT1-HeLaS3_INFS", 
                       "STAT1-HeLaS3_INFS-2", "ActivatedNucleosomes", 
                       "ActivatedNucleosomes.inout", "RestingNucleosomes", 
                       "RestingNucleosomes.inout", "Rest-CD4-HDAC1", 
                       "Rest-CD4-HDAC2", "Rest-CD4-HDAC3", "Rest-CD4-HDAC6", 
                       "Act-CD4-HDAC6", "Rest-CD4-CBP", "Rest-CD4-PCAF", 
                       "Rest-CD4-p300", "Rest-CD4-MOF", "Rest-CD4-Tip60", 
                       "Act-CD4-Tip60", "H3_3-HeLa", "p300-HeLa", 
                       "Brd4_2_1H_minus_100w_20s-Uwe", "Brd4_6_0H_minus_100w_20s-Uwe")
    intersect(display.order, c(meth.sets, act.sets, tf.sets))
}
