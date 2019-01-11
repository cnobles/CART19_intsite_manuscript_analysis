# Genomic Heatmap Maker for Integration Sites


# Script to run genomic heatmap maker for sites from DB

For list of samples and sample replicates generate heatmap
for a given reference genome and given integration site database:
```
Rscript genomic_heatmap_from_db.R sampleName_GTSP.csv -o heatmap --ref_genome hg18  --sites_group intSitesDev237
```

Group should be present in ~/.my.cnf.

File `sampleName_GTSP.csv` should have at least 2 columns: `sampleName` and
`GTSP`. `sampleName` is the name used in integration sites db. If column
`label` is given it will be used as a label in heatmap, otherwise GTSP column
will be used.  All samples with the same GTSP column merged together. See
geneTherapyPatientReportMaker 'check_patient_gtsp.R' script for instruction how
to generate `sampleName_GTSP.csv` for samples that are in speciment management
db.

## Database configuration file 

Configuration file should be in home directory and called .my.cnf,
(~/.my.cnf).

## Dependencies

intSiteRetriever
hiAnnotator
pipeUtils
colorspace
GCcontent

At present genomes hg18 and mm9 genomes added, if other genomes are
used add them to `genomicHeatmapMaker.R`.


## Implementation details: pipeUtils requirement

needs a column type which is hard-coded as "insertion" for integration site
and "match" for match random contol(mrc).

pipeUtils can only generate figures for 2 or more samples.


