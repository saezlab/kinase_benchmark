#!/usr/bin/env Rscript

# This path needs to be set correctly in order
# for the script to fins all necessary resources
install.location <- "INKA"

mapping.file <- paste0(install.location,"/valid_mapping.txt")
manning.kinase.file <- paste0(install.location,"/Manning.Kinases_20Jan2017.Rdata")

evquant.path <- paste0(install.location,"/evquant")

write_normalization_table <- function(mq.design = "INKA/experimentalDesignTemplate.txt",
                                   mq.evidence = "INKA/evidence.txt",
                                   normf.median = 1e7,
                                   mq.version = "1.5",
                                   out.normf = "INKA/normfactors.txt"

                                   ) {
    c.rawfile <- "Raw.file"

#   ---- experimentalDesignTemplate.txt file is read in
    tab.design <- read.delim(mq.design, header = TRUE)

    experiments <- unique(tab.design[,"Experiment"])


    no.raw <- nrow(tab.design)

    ## calculating normalization factor

    cat("Calculating normalization factors...")

    if (mq.version == "1.4") {
        evidence <- read.delim(mq.evidence, header = TRUE,
                       colClasses=c("Reverse" = "character",
                                    "Contaminant" = "character", "Experiment" = "character"))
    } else {
        evidence <- read.delim(mq.evidence, header = TRUE,
                                colClasses=c("Reverse" = "character",
                                             "Potential.contaminant" = "character",
                                             "Experiment" = "character"))
    }

    ## check for missing Experiment column
    exp.col <- "Experiment" %in% colnames(evidence)

    if (!exp.col) {
        cat("\nWhy isn't there a Experiment column in the evidence file?\n")
    } else {
        if (mq.version == "1.4") {
            evidence.ok <- subset(evidence, Reverse == "" &
                              Contaminant == "" & Intensity > 0 & MS.MS.Count > 0)
        } else {
            old.colnames <- colnames(evidence)
            old.colnames[grep("MS.MS..ount",colnames(evidence))] <- "MS.MS.Count"
            colnames(evidence) <- old.colnames
            evidence.ok <- subset(evidence, Reverse == "" & Potential.contaminant == "" &
                                    Intensity > 0 & MS.MS.Count > 0)
        }

        f <- rep(NA, length(experiments))

        for (i in 1:length(experiments)) {
            ex <- as.character(experiments[i])
            tmp <- subset(evidence.ok, Experiment == ex)
            if (nrow(tmp)>0) {
                f[i] <- median(tmp[,"Intensity"])
            }
        }

        if (normf.median > 0) {
            f.median <- normf.median
        } else {
            f.median <-  median(f, na.rm = TRUE)
        }

        normf <- f.median / f
        normf[is.na(normf)] <- 0

        write.table(cbind(experiment = gsub("[\\+-/ ()]", ".", as.character(experiments)), normf),
                file =  out.normf,
                sep='\t', row.names = F, quote = FALSE, qmethod = "double")

        cat(" Done.\n")
    }


}

#
# ------------
#

extract_mq_version <- function ( pp.file ){

    first.line <- readLines(pp.file, n=1)
    c.names <- unlist(strsplit(first.line,split="\t"))
    found <- grep("ontaminant",c.names)
    if (length(found)==0){
        version <- "1.5"
    }
    else {
        if (c.names[found] == "Potential contaminant" | c.names[found] == "Potential.contaminant"){
            version <- "1.5"
        }
        else {
            version <- "1.4"
        }
    }
    cat("Max Quant version detected as :",version,".\n")
    return(version)
}


#
# --------
#

to_kinase <- function(g.names_HGNC, kinases, sorted = TRUE) {

    tmp <- g.names_HGNC

    for (i in 1:length(tmp)) {
        if (tmp[i] != "") {
            gs <- unlist(strsplit(tmp[i], ";"))

            ind <- match(gs, kinases);

            a  <- NULL
            for (g in gs) {
                if (g %in% kinases) {
                    a <- c(a, g)
                }
            }

            if (is.null(a)) {
                tmp[i]  <-  ""
            }
            else {
                if (sorted) {
                    tmp[i]  <-  paste(sort(a), collapse = ";")
                } else {
                    tmp[i]  <-  paste(a, collapse = ";")
                }
            }
        }
    }

    return(tmp)
}

#
# --------
#

mapping_gn <- function(mapping, tab.pg.selected) {

    g.names <- NULL

    if (is.na(mapping)) {

        g.names <- as.character(tab.pg.selected[, 1]) # default gene names

    } else {

        mapping_file  <- read.delim(mapping, header = TRUE)
        rownames(mapping_file) <- mapping_file[, 1]

        g.names <- as.character(tab.pg.selected[, 2]) # protein ids

	for (i in 1:length(g.names)) {

            str <- as.character(g.names[i]);

            gn_tmp  <-  NULL;

            a <- unlist(strsplit(str,";"))

            if (length(a)>0) {
                for (j in 1:length(a)) {

                id <- a[j]

                    if (!is.na(id) && !is.null(mapping_file[id, 2]) && !is.na(mapping_file[id, 2])) {
                        if (as.character(mapping_file[id, 2]) != "") {
                            gn_tmp <- c(gn_tmp, as.character(mapping_file[id, 2]))
                        }
                    }
                }
            }
            if (is.null(gn_tmp)) {
                g.names[i]  <- ""
            } else {
                tmp2 <- NULL
                for (j in 1:length(gn_tmp)) {
                    if (is.na(match(gn_tmp[j], tmp2))) {
                        tmp2 <- c(tmp2, gn_tmp[j])
                    }
                }

                g.names[i] <- paste(tmp2, collapse=";")
            }
	}

    }

    return(g.names)
}

#
# --------
#

mapping_gn_HGNC <- function(mapping, tab.pg.selected) {

    g.names <- NULL

    if (is.na(mapping)) {

        g.names <- as.character(tab.pg.selected[, 1]) # default gene names

    } else {

        mapping_file  <- read.delim(mapping, header = TRUE)
        rownames(mapping_file) <- mapping_file[, 1]

        g.names <- as.character(tab.pg.selected[, 2]) # protein ids

	for (i in 1:length(g.names)) {

            str <- as.character(g.names[i]);

            gn_tmp  <-  NULL;

            a <- unlist(strsplit(str,";"))

            if (length(a)>0) {
                for (j in 1:length(a)) {

                    #id <- gsub("-.*", "", a[j]) # remove everything after and include "-""
                    id  <- a[j]

                    if (!is.na(id) && !is.null(mapping_file[id, 3]) && !is.na(mapping_file[id, 3])) {
                        if (as.character(mapping_file[id, 3]) != "") {
                            gn_tmp <- c(gn_tmp, as.character(mapping_file[id, 3]))
                        }
                    }
                }
            }
            if (is.null(gn_tmp)) {
                g.names[i]  <- ""
            } else {
                tmp2 <- NULL
                for (j in 1:length(gn_tmp)) {
                    if (is.na(match(gn_tmp[j], tmp2))) {
                        tmp2 <- c(tmp2, gn_tmp[j])
                    }
                }

                g.names[i] <- paste(tmp2, collapse=";")
            }
	}

    }

    return(g.names)
}

#
# -----
#

loadf <- function(f) {

        env <- new.env()
        obj <- load(f, env)
        env[[obj]]
}

#
# -----
#

write_pp_report <- function(file.in = "INKA/modificationSpecificPeptides.txt",
                            mq.design = "INKA/experimentalDesignTemplate.txt",
                            mq.evidence = "INKA/evidence.txt",
                            out.normf = "INKA/normfactors.txt",
                            mq.version = "1.5",
                            evquant.command = evquant.path,
                            kinase.file = "INKA/Manning.Kinases_12_2015.txt",
                            mapping = NA) {

    cat("Preparing pp-Report file.\n")

    if (mq.version == "1.4") {
        tab.pg <- read.delim(file.in, header = TRUE, colClasses=c("Reverse" = "character", "Contaminant" = "character"))

        tab.tmp0 <- subset(tab.pg, Reverse == "" & Contaminant == "" & Intensity > 0)

        tab.tmp <- subset(tab.pg, Reverse == "" & Contaminant == "" & Intensity > 0 & Phospho..STY. > 0)

    } else {

        tab.pg <- read.delim(file.in, header = TRUE, colClasses=c("Reverse" = "character", "Potential.contaminant" = "character"))

        tab.tmp0 <- subset(tab.pg, Reverse == "" & Potential.contaminant == "" & Intensity > 0)

        tab.tmp <- subset(tab.pg, Reverse == "" & Potential.contaminant == "" & Intensity > 0 & Phospho..STY. > 0)

    }


    cat("Mapping genes...")

    kinases <- as.character(loadf(kinase.file)[,1])

    tab.pg.selected <- tab.tmp[order(-tab.tmp[,"Intensity"]),]

    g.names <- paste0("=\"", mapping_gn(mapping, tab.pg.selected[, c("Gene.Names", "Proteins")]),  "\"")

    g.names_HGNC <- mapping_gn_HGNC(mapping, tab.pg.selected[, c("Gene.Names", "Proteins")])

    g.names_kinase <- to_kinase(g.names_HGNC, kinases)


    tab.design <- read.delim(mq.design, header = TRUE)

    experiments <- gsub("[\\+-/ ()]", ".", unique(tab.design[,"Experiment"]))

    col.int <- paste0("Intensity.", experiments)

    # normalize
    normf <- read.delim(out.normf, header = TRUE)
    rownames(normf) <- paste0("Intensity.", normf[,"experiment"])

    expr <- tab.pg.selected[, col.int]
    for (ex in colnames(expr)) {
        expr[, ex] <- expr[, ex] * normf[ex, "normf"]
    }
    colnames(expr) <- paste0("Norm.", colnames(expr))

    Protein.Groups <- paste0("=\"", as.character(tab.pg.selected[,"Protein.Groups"]), "\"")

    Gene.names <- paste0("=\"", as.character(tab.pg.selected[,"Gene.Names"]), "\"")

    Proteins <- tab.pg.selected[, c("Proteins")]

    cat(" Ready.\n")

    if (TRUE) {  ## create phosphopeptide MS/MS count

        cat("Creating phospho peptide count for peptide table ...\n")

        opl.pp.count <- tab.pg.selected[, col.int] * 0

        colnames(opl.pp.count) <- substr(colnames(opl.pp.count), 11, 1000)
# -----
# Need this one to restore column order in external
# C program output table.
        saved.column.order <- colnames(opl.pp.count)
#
        saved.column.order <- sub("^\\.+","X.",saved.column.order)
#	colnames(opl.pp.count) <- saved.column.order
# -----

        cat(colnames(opl.pp.count))
        cat("\n")

        tab.evidence <- read.delim(mq.evidence, header = TRUE)

        rownames(tab.evidence) <- tab.evidence[,"id"]

        old.v <- 0

# ------------------- ACCELERATED PART ---------------------------------------------
#
#
        write.table(tab.pg.selected,file="INKA/pp-input.txt",sep="\t",row.names=FALSE,quote=FALSE)
        run.command <- paste0(evquant.command," -o INKA/pp-output.txt -e ",mq.evidence," INKA/pp-input.txt")
        system(run.command)
#        file.remove("INKA/pp-input.txt")
        opl.pp.count <- read.delim("pp-output.txt",header=TRUE,sep="\t")
#        file.remove("INKA/pp-output.txt")
cat("Saved column order :",saved.column.order,"\n")
        opl.pp.count <- opl.pp.count[,saved.column.order]

# ----------------------------------------------------------------------------------

        colnames(opl.pp.count) <- paste0("OPL.ppCount.", colnames(opl.pp.count))

        expr <- cbind(expr, opl.pp.count)

        cat("Done.\n")
    }

    tab.out <- cbind(tab.pg.selected[, c("id",
                                         "Sequence",
                                         "Modifications",
                                         "Mass")],
                     Proteins,
                     list(Gene.names = g.names,
                          HGNCplus.gene.names = paste0("=\"", g.names_HGNC, "\""),
                          Sorted.kinase.names = paste0("=\"", g.names_kinase,  "\""),
                          MQ.gene.names = Gene.names),
                     tab.pg.selected[, c("Protein.Names",
                                         "Phospho..STY.",
                                         "MS.MS.Count",
                                         "Score",
                                         "Phospho..STY..site.IDs",
                                         col.int)],
                     expr)

    colnames(tab.out)[1] <- "ppModPeptide.ID"

    write.table(tab.out,
                file = "pp-Peptide-report.txt",
                sep = "\t",
                row.names = FALSE,
                quote = FALSE,
                qmethod = "double")

}

#
# ---------------
#

cleanup_files <- function (flist) {

    result <- file.remove(flist)
    if ( sum(result) != length(flist)){
        cat("Failed removing ",length(flist) - sum(result)," files.\n")
    }

}

#
# ---------------
#

mqver <-  extract_mq_version("INKA/modificationSpecificPeptides.txt")

write_normalization_table(mq.version=mqver)

write_pp_report(kinase.file = manning.kinase.file, mapping = mapping.file ,mq.version=mqver)

cat("CLeaning up temporary files.\n")
temp.files <- c("normfactors.txt" ,"pp-input.txt" ,"pp-output.txt")
cleanup_files(temp.files)
cat("Ready.\n")
