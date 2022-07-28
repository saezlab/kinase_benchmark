#!/usr/bin/env Rscript
#
# This script is a modified version for general 
# application of the inka pipeline to TMT encoded data.
#
# Author  : Alex Henneman
# Date    : Mon Jan 28 16:33:36 CET 2019
# Version : 1.0
#
# Usage:
# ======
# 
# 1) Prepare (header containg) tab separated table with TMT
#    postfix (What's appended to "Reporter.intensity.corrected." in the
#    modificationSpecificPeptides.txt table) and the desired sample name.
# 2) Run in MaxQuant search results directory
# 3) Run inka.R script
# ------


cat("Reading modpep file...")
if(!exists("mp")){
    mp <- read.delim("modificationSpecificPeptides.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
}
cat(" Ready.\n")

# This is the format of this file
# > colnames(rep.labels)
# [1] "TMT.Tag"    "Sample.Tag"
rep.labels <- read.table("sample_reporter_table.txt",header=TRUE,colClasses=c("character","character"),stringsAsFactors=FALSE)
cat("Processing ",nrow(rep.labels)," TMT channnels.\n")

tmt.sample.tag <- "Reporter.intensity.corrected."
inka.sample.tag <- "OPL.ppCount."
# Here is where the quantifications are: Reporter.intensity.corrected.XXX
for (i in seq_along(rep.labels$TMT.Tag)){

    tmt.tag    <- rep.labels$TMT.Tag[i]
    sample.tag <- rep.labels$Sample.Tag[i]

    old.sample.column <- paste0(tmt.sample.tag,tmt.tag)
    new.sample.column <- paste0(inka.sample.tag,sample.tag)
    cat("Converting ",old.sample.column," to ",new.sample.column,"\n")

    mp[[ new.sample.column ]] <- mp[[ old.sample.column ]]
}

## Going to filter out input data
## get rid of these :"Reverse"                   "Potential.contaminant"
row.select <- mp$Reverse != "+" & mp$Potential.contaminant != "+" 
cat("Non reverse contaminant rows: ",sum(row.select)/nrow(mp)*100,"%\n")
mp <- mp[row.select,]
# Get rid of zero peptide rows
inka.cols <- grep(inka.sample.tag,colnames(mp))
row.nz <- rowSums(mp[,inka.cols])>0
cat("Non-zero peptide rows: ",sum(row.nz)/nrow(mp)*100,"%\n")
mp <- mp[row.nz,]
#
outfile.name <- "pp-Peptide-report.txt"
mp$ppModPeptide.ID <- mp$id
#
cat("Writing output to file ",outfile.name,"\n")
write.table(mp,file=outfile.name,quote=FALSE,row.names=FALSE,sep="\t")
