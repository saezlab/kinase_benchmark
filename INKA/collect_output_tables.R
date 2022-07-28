#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

generate.tables.version <- "0.0.4"
cat("This generate_inka_tables.R version :",generate.tables.version,"\n")

prefix.tag <- ""

if (length(args)>0){
    prefix.tag <- args[1]
    cat("Using prefix tag:",prefix.tag,"\n")
}

mqread <- function ( fname ){
	read.delim(fname,header=TRUE,stringsAsFactors=FALSE)
}

translate_sample_tags <- function (sample.tags,label.file){
	
	labs <- read.delim(label.file,header=FALSE,sep='\t',stringsAsFactors=FALSE)
    labs[[ncol(labs)]] <- make.names(labs[[ncol(labs)]])
	last.col <- ncol(labs)
	for ( i in seq_along(sample.tags)){
		old.tag <- sample.tags[i]
		old.row.index <- match(old.tag,labs[[1]])
        if (!is.na(old.row.index)){
		    new.tag <- labs[[last.col]][old.row.index]
		    sample.tags[i] <-  sub("\\.+$","",make.names(new.tag))
		    cat(" * Replacing name ",old.tag," in row ",old.row.index," by ",new.tag,"\n")
        } else {
            cat(" * Tag ",old.tag," not found, therefore left unchanged\n")
        }
	}
	return(sample.tags)
}

sample.name.source <- "pp-Peptide-report.txt"
pp <- mqread(sample.name.source)
all.col.names <- colnames(pp)
common.sample.prefix <- "OPL.ppCount."
sample.tags <- sub(common.sample.prefix,"",all.col.names[grep(common.sample.prefix,all.col.names)])
cat("Original names: ",sample.tags,"\n")

label.file <- "labels.txt"
if (file.exists(label.file)){
    sample.tags <- translate_sample_tags(sample.tags,label.file)
} 

raw <- list(length(sample.tags))

cat("Using names: ",sample.tags,"\n")

for (i in seq_along(sample.tags)){
    cat("Importing sample data from ",sample.tags[i],"\n")
    input.filename <- paste0("inka_save_",sample.tags[i],".txt")
    raw[[i]] <- mqread(input.filename)
#file.remove(input.filename)
}

names(raw) <- sample.tags

cat("Gathering kinase names.\n")
all.kinase.names <- NULL
for (i in seq_along(sample.tags)){
    all.kinase.names <- unique(c(all.kinase.names,raw[[i]]$Kinase))
}

first.part <- data.frame(Kinase=all.kinase.names,stringsAsFactors=FALSE)
n.kins <- length(all.kinase.names)

kin.all  <- first.part
act.all  <- first.part
psp.all  <- first.part
nwk.all  <- first.part
inka.all <- first.part

for (i in seq_along(sample.tags)){
    cat("Appending sample ",sample.tags[i],"\n")
    df <- raw[[i]]
    kin.all <- merge(kin.all,data.frame(Kinase=df$Kinase,value=df$Kinome,stringsAsFactors=FALSE),by="Kinase")
    tmp <- colnames(kin.all)
    tmp[length(tmp)] <- paste0("Kinome_",sample.tags[i])
    colnames(kin.all) <- tmp
#
    act.all <- merge(act.all,data.frame(Kinase=df$Kinase,value=df$ActLoop,stringsAsFactors=FALSE),by="Kinase")
    tmp <- colnames(act.all)
    tmp[length(tmp)] <- paste0("ActLoop_",sample.tags[i])
    colnames(act.all) <- tmp
#
    psp.all <- merge(psp.all,data.frame(Kinase=df$Kinase,value=df$PSP,stringsAsFactors=FALSE),by="Kinase")
    tmp <- colnames(psp.all)
    tmp[length(tmp)] <- paste0("PSP_",sample.tags[i])
    colnames(psp.all) <- tmp
#
    nwk.all <- merge(nwk.all,data.frame(Kinase=df$Kinase,value=df$NWK,stringsAsFactors=FALSE),by="Kinase")
    tmp <- colnames(nwk.all)
    tmp[length(tmp)] <- paste0("NWK_",sample.tags[i])
    colnames(nwk.all) <- tmp
#
    inka.all <- merge(inka.all,data.frame(Kinase=df$Kinase,value=df$Score,stringsAsFactors=FALSE),by="Kinase")
    tmp <- colnames(inka.all)
    tmp[length(tmp)] <- paste0("InKA_",sample.tags[i])
    colnames(inka.all) <- tmp
}

# ---------------------------------------------------

select.above.threshold <- function ( tab, thresh ){
    selected.rows <- rowSums(tab[,-1] >= thresh) > 0
    tab <- tab[selected.rows,]
}

# ---------------------------------------------------

remove.zero.kinase.rows <- function (tab){

    non.zero.rows <- rowSums(tab[,-1]) > 0
    tab <- tab[non.zero.rows,]
    return(tab)
}

# ---------------------------------------------------

kin.all  <- remove.zero.kinase.rows(kin.all)
act.all  <- remove.zero.kinase.rows(act.all)
psp.all  <- remove.zero.kinase.rows(psp.all)
nwk.all  <- remove.zero.kinase.rows(nwk.all)
inka.all <- remove.zero.kinase.rows(inka.all)

# 
# ---------------------------------------------------
#

save.sample.table <- function ( tab,name ){

    write.table(tab,file=name,sep="\t",row.names=FALSE)
    cat("Output written to file ",name,".\n")

}

# 
# ---------------------------------------------------
#

cat("Saving tables...\n\n")
kin.save.file  <- paste0(prefix.tag,"kinome_table.tsv")
act.save.file  <- paste0(prefix.tag,"actloop_table.tsv")
psp.save.file  <- paste0(prefix.tag,"psp_table.tsv")
nwk.save.file  <- paste0(prefix.tag,"nwk_table.tsv")
inka.save.file <- paste0(prefix.tag,"inka_table.tsv")

save.sample.table(kin.all,kin.save.file)
save.sample.table(act.all,act.save.file)
save.sample.table(psp.all,psp.save.file)
save.sample.table(nwk.all,nwk.save.file)
save.sample.table(inka.all,inka.save.file)
cat("\nReady.\n")

