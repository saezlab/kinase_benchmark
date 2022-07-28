#!/usr/bin/env Rscript

cat("Converting pp-file to intensities.\n")

pp.filename <- "pp-Peptide-report.txt"
outfile.name <- "pp-Peptide-report_intensities.txt"

pp <- mqread(pp.filename)

count.col.tag     <- "OPL.ppCount."
intensity.col.tag <- "Norm.Intensity."

count.col.names <- colnames(pp)[grep(count.col.tag,colnames(pp))]
sample.names <- sub(count.col.tag,"", count.col.names)

for (i in seq_along(sample.names)){
    src.name  <- paste0(intensity.col.tag,sample.names[i] )
    dest.name <- paste0(count.col.tag,sample.names[i] )
    cat(" * Going to copy intensitie values from ",src.name," to ",dest.name,"\n")
    pp[[ dest.name ]] <- pp[[ src.name ]]
}


write.table(pp,file = outfile.name,sep = "\t",row.names = FALSE,quote = FALSE,qmethod = "double")
cat("Written results to file ",outfile.name,"\n")
