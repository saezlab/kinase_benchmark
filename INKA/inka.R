#!/usr/bin/env Rscript 

# Set this variable to absolute path of
# installation directory
general.resources.dir <- "."

inka.script.version <- "0.5.1"

args = commandArgs(trailingOnly=TRUE)

# These are the default values
network.link.threshold <- ">2"
inka.pdf.output.filename <- "inka_out.pdf"
tiox.mode <- FALSE

if (length(args)>0){
    match.tiox <- match("--tiox",args)
    if (!is.na(match.tiox)){
        cat("Running in TiOx mode.\n")
        network.link.threshold <- ">10"
        args <- args[ - match.tiox ]
        tiox.mode <- TRUE
    }
#
#   These should be the left over 
#   flags arguments
    if (length(args)>0){
        inka.pdf.output.filename <- args[1]
    }

} 
cat("Network link thershold equal to :",network.link.threshold,"\n")


inkascore.version <- "0.2.1"

pp.filepath <- "./pp-Peptide-report.txt"
ps.filepath <- "./Phospho (STY)Sites.txt"

require(tools)
library(data.table)
library(splitstackshape)
library(stringr)
library(network)
library(gplots)

options(scipen=10)

loadf <- function(f) {

        env <- new.env()
        obj <- load(f, env)
        env[[obj]]
}
        

#
# -----------------------------------------------------------------------------
# 

relabel <- function(x, snames,data.type) {
                 pref <- paste0(data.type, ".")
                 subl <- c("_Kinome", "_ActivationLoop", 
                           "_PhosphoSitePlus", "_NetworKIN")
                 names(x) <- paste0(pref, snames, rep(subl, each = length(x)/4))
                 return(x)
}
    

#
# -----------------------------------------------------------------------------
# 

getData <- function(TABLE, dtype ) {

        dt         <- copy(TABLE)
        kpcol      <- match(c("KINASE", "Phospho.Mods"), names(dt))
        scol       <- grep(dtype, names(dt))
        numSamples <- length(scol)
        samples    <- names(dt) [scol]
        PLOT_DATA  <- list()

        for (item in 1:numSamples) {
                
            if(item > 1) dt <- copy(TABLE)
            sel <- c(kpcol[1], scol[item], kpcol[2])
            dt  <- dt [, .SD, .SDcols = sel]
            setnames(dt, c("KINASE", "SIGNAL", "MODS"))
            dt [, SIGNAL := SIGNAL * MODS]
            dt [, MODS   := NULL]
            PLOT_DATA[[item]] <- dt
        }
        names(PLOT_DATA) <- samples

        return(PLOT_DATA)
}
        

#
# -----------------------------------------------------------------------------
# 
        
getData2 <- function(TABLE, dtype ) {
        
        dt.orig <- copy(TABLE)
        dt.orig <- dt.orig [! duplicated(paste(Ppeptide.Id, Position.In.Peptide, KINASE))]
        dt.orig [, N.Sites := .N, by = .(Ppeptide.Id, KINASE, SUBSTRATE)]

        kpcol      <- match(c("KINASE", "Phospho.Mods", "N.Sites"), names(dt.orig))
        scol       <- grep(dtype, names(dt.orig))
        numSamples <- length(scol)
        samples    <- names(dt.orig) [scol]
        PLOT_DATA  <- list()

        for (item in 1:numSamples) {
                
            dt <- copy(dt.orig)

            sel <- c(kpcol[1], scol[item], kpcol[2:3])
            dt <- dt [, .SD, .SDcols = sel]
            setnames(dt, c("KINASE", "SIGNAL", "MODS", "N.SITES"))
            
            dt [, SIGNAL := SIGNAL * MODS/N.SITES ]
            dt [, `:=` (MODS = NULL, N.SITES = NULL)]
            PLOT_DATA[[item]] <- dt
        }

        names(PLOT_DATA) <- samples

        return(PLOT_DATA)
}

#
# -----------------------------------------------------------------------------
# 
       
       
foldSort <- function(x, folds = 4) { 
                
        xlen <- length(x)
        flen <- xlen/folds
        
        if(flen %% 1 != 0) {
            stop("Nummber of sort items is not a multiple of the same number.")
        }

        xsrt <- as.integer(t(matrix(1:xlen, nrow = flen)))
                       
        return( x [xsrt] )
}

#
# -----------------------------------------------------------------------------
# 

        
process_P <- function(report_PP) {                                               

    tmp <- names(report_PP)
    tmp <- sub("\\.+$", "", gsub(" +|\\.+", ".", tmp))
    tmp [tmp == "ppModPeptide.ID"]      <- "Ppeptide.Id"
    tmp [tmp == "Phospho.STY"]          <- "Phospho.Mods"
    tmp [tmp == "Phospho.STY.site.IDs"] <- "Psite.Id"
    setnames(report_PP, tmp)
	
    old  <- match("Psite.Id", names(report_PP))
    last <- ncol(report_PP)
    srt  <- c(1, old, 2:(old-1), (old+1):last)
	
    report_PP <- report_PP [, .SD, .SDcols = srt]
    report_PP [, Gene := Gene.Names]

    report_PP$Ppeptide.Id <- as.character(report_PP$Ppeptide.Id)
    
    return(report_PP)
}

#
# -----------------------------------------------------------------------------
# 

isosort <- function(accs) {
            
    accs      <- unlist(strsplit(accs, split = ";")) 
    tmp       <- accs  
    sel       <- ! grepl("-", tmp)
    tmp [sel] <- paste0(tmp [sel], "-0")   
    ord       <- order(
                 as.integer(sapply(strsplit(tmp, "-"),"[[", 2)))
    accs      <- paste(accs[ord], collapse = ";")                                 
    return(accs)
} 

#
# -----------------------------------------------------------------------------
# 

process_S <- function(maxQ_PS,contaminant) {
	
    
#

suppressWarnings(          
    maxQ_PS [, Sorted.Accessions := sapply(Proteins, isosort)]
)   
    colpos <- match(c("id", "Mod. peptide IDs", "Sequence window", 
                      "Amino acid", "Localization prob", 
                      "Sorted.Accessions", "Proteins", 
                      "Positions within proteins", "Reverse", 
                       contaminant),                          
                    names(maxQ_PS))
    maxQ_PS <- maxQ_PS [,  .SD, .SDcols = colpos]
    setnames(maxQ_PS, 
             c("Psite.Id", "Ppeptide.Id", "SeqWindow.31AA", 
               "Amino.Acid", "Localization.probability", 
               "Sorted.Accessions", "Accession", 
               "Position.In.Accession", "Reverse.hit", 
               "Contaminant.hit"))
    setkey(maxQ_PS, Reverse.hit, Contaminant.hit)                             
    maxQ_PS <- maxQ_PS [.("","")]                                            
    maxQ_PS [, `:=` (
        Reverse.hit = NULL, 
        Contaminant.hit = NULL, 
        Spacer = "ANNOTATION >>",
        ClassI = ifelse(Localization.probability > 0.75, "+", "-"))]
    
    return(maxQ_PS)                                                           
} 

#
# -----------------------------------------------------------------------------
# 

accessionise <- function( classI_PS, UniProt_data ) {
        
    accWise_classI_PS <- cSplit(classI_PS,                                  
	                        splitCols    = c("Accession", 
					       "Position.In.Accession"), 
	                        sep          = ";", 
                                direction    = "long", 
			        type.convert = FALSE)
    accWise_classI_PS [, `:=`(
        SeqWindow.15AA  = substr(SeqWindow.31AA, 9, 23),
        Accession_AApos = paste0(Accession,"-",
                                 Amino.Acid, Position.In.Accession))]
    old               <- match("Accession", names(accWise_classI_PS))
    accWise_classI_PS <- merge(x  = accWise_classI_PS,                                      
                               y  = UniProt_data,                                                   
                               by = "Accession", all.x = TRUE)               
    last              <- ncol(accWise_classI_PS)
    srt               <- c(2:old, 1, (old+1):last)
    accWise_classI_PS <- accWise_classI_PS [,  .SD, .SDcols = srt]
    
    return(accWise_classI_PS)
}

#
# -----------------------------------------------------------------------------
# 

filtersort <- function ( accWise_classI_PS,PSP_human,PSP_KSR_human ) {
        
    XinY <- function(x,y) ifelse(x %in% y, "+", "-")
    alfN <- function(x) { y <- utf8ToInt(substr(x, 1, 1)) - 64
                          y <- ifelse(y == 16, 100, y) }
    indx <- function(x,y) match(x, unlist(strsplit(y, split=";")))
                
    f_accWise_classI_PS <- accWise_classI_PS [! is.na(UniProt.Gene)]
    
    f_accWise_classI_PS [, Gene.multiplicity := 
        length(unique(UniProt.Gene)), by = Psite.Id ]
    
    f_accWise_classI_PS [, `:=` (
        Site.in.PSP = 
            XinY(Accession_AApos, PSP_human$Accession_AApos),
        Window.in.PSP = 
            XinY(SeqWindow.15AA, PSP_human$SeqWindow.15AA),
        Site.in.PSP_KSR = 
            XinY(Accession_AApos, PSP_KSR_human$Accession_AApos),
        Window.in.PSP_KSR =
            XinY(SeqWindow.15AA, PSP_KSR_human$SeqWindow.15AA),
        Primary.Accession = 
            ifelse(Accession == UniProt.Entry, "+", "-"),
        Alphabet.Pos = 
            sapply(Accession, alfN),
        Accession.Idx = 
            mapply(indx, Accession, Sorted.Accessions)
        )]

    setorder(f_accWise_classI_PS, 
	                                Psite.Id, 
	                                UniProt.Gene, 
	                                Window.in.PSP,
	                                Window.in.PSP_KSR,
	                                Site.in.PSP, 
	                                Site.in.PSP_KSR,
	                               -UniProt.Score,
	                                Primary.Accession,
	                               -Alphabet.Pos,
	                                Accession.Idx)
    
    f_accWise_classI_PS [, `:=` (
        Alphabet.Pos  = NULL, 
        Accession.Idx = NULL)]
    
    return(f_accWise_classI_PS)
}

#
# -----------------------------------------------------------------------------
# 
       
posinpep <- function(win, pep) {
            
    pos <- NULL
    site.pos <- ceiling(nchar(win)/2)
    site.aa <- substr(win, site.pos, site.pos)
            
    matches <- unlist(
               gregexpr(paste0('(?=', site.aa, ')'), pep, perl=T))
                              
    for (cand in matches) {
                
        win_l <- substr(win, 1, site.pos)
        pep_l <- substr(pep, 1, cand)
        win_r <- substr(win, site.pos, nchar(win))
        pep_r <- substr(pep, cand, nchar(pep))
                
        nwl <- nchar(win_l)
        npl <- nchar(pep_l)
        nwr <- nchar(win_r)
        npr <- nchar(pep_r)
                
        if(npl>nwl) {
                pep_l <- substr(pep_l, npl-nwl+1, npl)
        } else {
                win_l <- substr(win_l, nwl-npl+1, nwl)
        }
                
        if(npr>nwr) {
                pep_r <- substr(pep_r, 1, nwr)
        } else {
                win_r <- substr(win_r, 1, npr)
        }
                
        if(pep_l != win_l || pep_r != win_r) {
                next 
        } else {
                pos <- c(pos, cand)
        }
    }
                
    if(is.null(pos)) pos <- "No match found!"
    
    return(paste(as.character(pos), collapse = ";"))
}

#
# -----------------------------------------------------------------------------
# 

pepulate <- function ( f_accWise_classI_PS,report_PP ) {                                    
                        
#   NOTE That here we only keep the first entry we come across
    nonRed_classI_PS <- f_accWise_classI_PS [
                           !duplicated(paste(Psite.Id, toupper(UniProt.Gene)))]  
#   Here we filter only those Ppeptide.IDs that are also in our
#   report_PP$Ppeptide.Id
    nonRed_classI_PS_PP <- cSplit(nonRed_classI_PS, splitCols = "Ppeptide.Id", 
				                  sep       = ";", direction = "long")
#    NOTE needed for driver script
    nonRed_classI_PS_PP$Ppeptide.Id <- as.character(nonRed_classI_PS_PP$Ppeptide.Id)
    nonRed_classI_PS_PP <- nonRed_classI_PS_PP [
                               Ppeptide.Id %in% report_PP$Ppeptide.Id]
#   Here we merge nonRed_classI_PS_PP with a select set of columns
#   from report_PP by Ppeptide.Id
    cols <- grep( "Ppeptide.Id|Sequence|Phospho.Mods|Spect", names(report_PP))
    nonRed_classI_PS_PP <- merge(x  = nonRed_classI_PS_PP, 
                                 y  = report_PP [, .SD, .SDcols = cols], 
                                 by = "Ppeptide.Id")
#   Add a column Position.In.Peptide to merge and reorder
#   columns
    nonRed_classI_PS_PP [, Position.In.Peptide := 
        mapply(posinpep, SeqWindow.15AA, Sequence)]              
    behind <- match("Phospho.Mods", names(nonRed_classI_PS_PP))
    old    <- ncol(nonRed_classI_PS_PP)
    srt    <- c(1:behind, old, (behind+1):(old-1))

    nonRed_classI_PS_PP <- nonRed_classI_PS_PP [, .SD, .SDcols = srt]
    
    return(nonRed_classI_PS_PP)

}

#
# -----------------------------------------------------------------------------
# 

determine_contaminant_tag <- function (maxQ){

    c.tag <- "Contaminant"
    if (maxQ != "1.4") c.tag <- "Potential contaminant"

    return(c.tag)
}

#
# -----------------------------------------------------------------------------
# 
    
read_mq_ps_file <- function ( path,maxQ){

    contaminant <- determine_contaminant_tag(maxQ)

    selected.columns <- c("Proteins", "Positions within proteins", 
                          "Localization prob", "Amino acid", 
                          "Sequence window", "Reverse", 
                          contaminant,"id", "Mod. peptide IDs")
    ps.full <- fread(path)
    ps <- ps.full[,selected.columns, with=FALSE]
#    This does not seem to work on some files.
#    ps <- fread(path, select = c("Proteins", "Positions within proteins", 
#                                            "Localization prob", "Amino acid", 
#                                            "Sequence window", "Reverse", 
#                                             contaminant,
#                                            "id", "Mod. peptide IDs"))
    return (ps)
}

# -----------------------------------------------------------------------------

filter_annotate <- function( read.ppeptides,read.psites,UniProt_data,PSP_human,PSP_KSR_human,maxQ){


        contaminant  <- determine_contaminant_tag(maxQ)
        
	    report_PP           <- copy(read.ppeptides)
	
        report_PP           <- process_P(report_PP)

        geneWise_report_PP  <- cSplit(report_PP,                                 
		        	      splitCols = "Gene", 
			                    sep = ";", 
			              direction = "long",
		                   type.convert =  FALSE)
	
        report_PP [, Gene := NULL]
        
#       independent of above

        maxQ_PS             <- copy(read.psites)
        
        noRev_noCon_PS      <- process_S(maxQ_PS,contaminant)                                     
        
        setkey(noRev_noCon_PS, ClassI)
        
        classI_PS           <- noRev_noCon_PS ["+"]

        accWise_classI_PS   <- accessionise(classI_PS,UniProt_data)                           
        
        f_accWise_classI_PS <- filtersort(accWise_classI_PS,PSP_human,PSP_KSR_human)
         
        nonRed_classI_PS_PP <- pepulate(f_accWise_classI_PS,report_PP)
   
        cat("Finished data filtering & annotation.\n")
        
        retval <- list(report_PP,geneWise_report_PP,nonRed_classI_PS_PP)
        names(retval) <- c("report_PP","geneWise_report_PP","nonRed_classI_PS_PP")
        return(retval)
}

#
# -----------------------------------------------------------------------------
# 


tabulate_kinome <- function (PP,Kinases,mapgen) {                                  

	        gscolP <- match(mapgen, names(PP))
	        gscolK <- match(mapgen, names(Kinases))

	        names(PP) [gscolP] <- "mappedGS"
    PP <- PP [Gene == mappedGS | ! mappedGS %in% Gene]
    PP <- PP [! duplicated(paste(Ppeptide.Id, mappedGS))]
    
    kinase.GS <- Kinases [[gscolK]]

    KINOME_TABLE <- PP [mappedGS %in% kinase.GS]
    KINOME_TABLE [, KINASE := mappedGS]                                   
    KINOME_TABLE [, KIN_SEQ := paste(KINASE, Sequence, sep="_")]
    last <- ncol(KINOME_TABLE)
    srt  <- c(last-1:0, 1:(last-2))
	        KINOME_TABLE <- KINOME_TABLE [, .SD, .SDcols = srt] 
	    
    return(KINOME_TABLE)

}

#
# -----------------------------------------------------------------------------
# 

# NOTE that PP here is a read-in genewise thingy
generate_kinome_pdata <- function(PP,Kinases,mapgen){
        
    KINOME_TABLE <- tabulate_kinome(PP,Kinases,mapgen)

	return(KINOME_TABLE)                     
}

#
# -----------------------------------------------------------------------------
# 


tabulate_act_loop <- function (PP,Kinases,KINOME_TABLE,mapgen) {                                          

    setnames(PP, gsub(" ", ".", names(PP)))
    PP <- PP [! is.na(Loop.Sites)]
    PP [, KINASE := sapply(strsplit(Gene.Name, split=" "), "[[",1)]
    
	        gscolK    <- match(mapgen, names(Kinases))
	        kinase.GS <- Kinases [[gscolK]]
    notFound  <- PP [! KINASE %in% kinase.GS]

    if(nrow(notFound) != 0) {
        
        print(data.frame(notFound))
        stop("unmapped kinases\n\n",
             "KINASE entries in the subtable above were not in the ",
             "'Kinases' table.\nPlease edit the pertinent kinase gene ", 
             "symbols ('Gene Name') in the 'phomics_results'\ntable to", 
             " ones that are among the HGNC-mapped gene symbols in the", 
             " 'Kinases' table.")
    }

    PP [, KIN_SEQ := paste(KINASE, toupper(Peptide), sep = "_")]
    loop.codes     <- PP [, KIN_SEQ]
    ACT_LOOP_TABLE <- KINOME_TABLE [KIN_SEQ %in% loop.codes]

    return(ACT_LOOP_TABLE)
}

#
# -----------------------------------------------------------------------------
# 

generate_act_loop_pdata <- function(PP,KINOME_TABLE,Kinases,mapgen,dtype) {
        
        cat("Starting activation loop plot data generation ...")
        
        ACT_LOOP_TABLE <- tabulate_act_loop(PP,Kinases,KINOME_TABLE,mapgen)   
        
        cat("Ready.\n")
        
	return(ACT_LOOP_TABLE)
}

#
# ----------------------------------
#

generate_PSP_pdata <- function( nonRed_classI_PS_PP,PSP_KSR_human, mapgen ,dtype) {
        
    cat("Starting PSP plot data generation...")

    PSP_PS <- merge(x  = PSP_KSR_human, 
                    y  = nonRed_classI_PS_PP,                                           
                    by = "Accession_AApos")
    
    cols <- c(mapgen, "UniProt.Gene", "Ppeptide.Id", "Phospho.Mods",
              "Sequence", "Position.In.Peptide", "Psite.Id", 
              "SeqWindow.15AA.y", "Accession_AApos", dtype)
    srt  <- unlist(sapply(cols, function(x) grep(x, names(PSP_PS))))
    
    PSP_PS_TABLE <- PSP_PS [, .SD, .SDcols = srt]
    
    tmp <- c("KINASE", "SUBSTRATE", names(PSP_PS_TABLE)[-(1:2)])
    tmp <- sub("\\.y", "", tmp)
    setnames(PSP_PS_TABLE, tmp)
        cat(" Ready.\n")
        
	return(PSP_PS_TABLE)  
}

#
# ----------------------------------
#

tabulate_NWK <- function (TOP_NWK,filtData,dtype) {                   

    old <- match("Kinase.STRING.ID", names(TOP_NWK))
    names(TOP_NWK)[old] <- "Kinase.Ensembl"
    TOP_NWK <- TOP_NWK [, .SD, .SDcols = c(2:ncol(TOP_NWK), old)]
    #TOP_NWK <- TOP_NWK [, .SD, .SDcols = c(setdiff(1:ncol(TOP_NWK), old),old) ]

    NWK_PS  <- merge(x  = TOP_NWK, 
                     y  = filtData,                                       
                     by = "Accession_AApos")
    
    cols <- c("Kinase.HGNC", "UniProt.Gene", "Ppeptide.Id", 
              "Phospho.Mods", "Sequence", "Position.In.Peptide",
              "Psite.Id", "Accession_AApos", "SeqWindow.15AA", dtype)
    srt  <- unlist(sapply(cols, function(x) grep(x, names(NWK_PS))))
    
    NWK_PS_TABLE <- NWK_PS [, .SD, .SDcols = srt]
    
    temp <- c("KINASE", "SUBSTRATE", names(NWK_PS_TABLE)[-(1:2)])
    setnames(NWK_PS_TABLE, temp)

    return(NWK_PS_TABLE)                                                
}

#
# ----------------------------------
#

generate_NWK_pdata <- function(TOP_NWK,filtData,dtype,
                               mapgen){
        
        cat("Starting NWK plot data generation...")

        NWK_PS_TABLE <- tabulate_NWK(TOP_NWK,filtData,dtype) 
        
        cat(" Ready.\n")
        
	return(NWK_PS_TABLE)
}

#
# ----------------------------------
#

bplot <- function(DATA,name,kinase=NULL) {     
        
    font.scale <- 0.5
    dt <- copy(DATA)
        
    dt [, SIGNALsum := sum(SIGNAL), by = KINASE]
    setorder(dt, -SIGNALsum, KINASE)
    dt [, Rank := sapply(KINASE, function(x) match(x,unique(KINASE)))]
    TOP20       <- dt  [Rank<=20]
    MAXN        <- TOP20 [, .N, by = KINASE] [, max(N)]
    kinames     <- TOP20 [, KINASE]
    data.vector <- TOP20 [, SIGNAL]
    
    split.data  <- split(data.vector, kinames)

    p <- sapply(split.data, function(x) c(x, rep(0, MAXN - length(x))))
    S <- colSums(p)
    M <- max(S)
    p <- p [, order(S)]
    p <- apply(p, 2, sort, decreasing = TRUE)
    
    suppressWarnings(
    barplot(p, 
            space = 0, horiz = TRUE, yaxs = "i", axes = FALSE, 
            cex.names = 0.8*font.scale, line = -0.6, las = 1, col = "grey90")
    )
    axis(side = 1, at = floor(c(0, M)), labels = FALSE, tcl = 0.25)                
    axis(side = 1, padj=0,cex.axis = 0.8*font.scale)                
    
    mtext(name, side = 1, at = 0.6 * M, line= -1.5 , cex = font.scale, font = 2)

    if(! is.null(kinase)) {
            
        for (kin in kinase) {
                
            if(is.na(kin)) next
        
            p2        <- matrix(0, nrow(p), ncol(p))
            highlight <- colnames(p) %in% kin
            idx       <- which(kinase == kin)
        
            for (rw in 1:nrow(p2)) {
                
                p2 [rw, highlight] <- p [rw, highlight]
                
            }
            
            barplot(p2, 
                    space = 0, horiz = TRUE, yaxs = "i", 
                    axes = FALSE, cex.names = 0.8*font.scale, las = 1, 
                    col = cl[idx], add = TRUE)
        }
    }

    return(invisible(S))
}  

#
# ----------------------------------
#

kinase_condense_data <- function ( full.data ,quant.name=c("KINASE","KINOME","ACTLOOP","PSP","NWK")){

    dt <- copy(full.data)
    kins    <- unique(unlist(lapply(dt, function(x) x [, KINASE])))
    cat("Condensing data into ",length(kins)," kinases.\n")
    retval <- list(KINASE=kins)
    for (i in 1:4){
        tmp <- dt[[i]]
        tmp <- tmp[, sum(SIGNAL) ,by = KINASE]
        colnames(tmp) <- c("KINASE",quant.name[i+1])
        retval <- merge(retval,tmp,by="KINASE",all=TRUE)
    }
    retval[is.na(retval)] <- 0
    retval <- data.table(retval)

    return (retval)
}

#
# ----------------------------------
#

make_kinase_bar_plot <- function( DATA,name,frame ) {
        
        
    cat("Starting single sample kinase bar ploting ...")

    par(mar=c(1.5,1,1,1))
        
    width <- frame[2] - frame[1]
    singl.width <- width/4
    prefix <- c("Kinome","ActLoop","PSP","NWK")
    for (idx in 1:length(DATA)) { 
        x.start <- frame[1]+(idx-1)*singl.width
        x.stop  <- frame[1]+idx*singl.width
        new.frame <- c(x.start,x.stop,frame[3],frame[4])
        par(fig=new.frame)
        bplot(DATA [[idx]],paste0(name,"\n",prefix[idx]))
        par(new=TRUE)
    }

    cat("Ready.\n")
        
}


#
# ----------------------------------
#


make_inka_plot <- function ( pdt, sample.name, mapgen,Kinases,frame, 
                        tiox.flag=FALSE,
                        COLkinases     = NULL, BARCOLkinases  = NULL, clb            = NULL,
                        cl             = c("blue", "black", "grey60", "grey80")) {
        
        font.size.scale <- 0.5
        par(mar=c(2,1,1,1))

        cat("Starting INKA plotting ...with pdt ",nrow(pdt)," rows\n")

        frame.div <- frame[1] + 0.25*(frame[2]-frame[1])
        par(fig=c(frame[1],frame.div,frame[3],frame[4]),new=TRUE)
        
        barKins <- Kinases [Kinase_without_Substrate > 0, eval(as.name(mapgen))]
        cat(length(barKins)," kinases out of PSP/NWK\n")
        
        par(oma = c(0, 0.1, 0, 0.1)) 
        
        pdt [is.na(pdt)] <- 0               
        setnames(pdt, c("Kinase", "Kinome", "ActLoop", "PSP", "NWK"))
        
        kdt <- pdt [Kinase %in% barKins & Kinome >= 2]
        cat("Out of PSP/NWK kinases for this sample",nrow(kdt)," \n")
        kdt [, Kin := Kinome + ActLoop]
        setorder(kdt, -Kin)
        
        p <- kdt[20:1, {tmp <- Kin; names(tmp) <- Kinase; tmp}]
        p [is.na(p)] <- 0
        
        suppressWarnings(
        barplot(p, 
                horiz = TRUE, space = 0, yaxs = "i", axes = FALSE,
                las = 1, line = -0.5, col = "grey95", cex.names = 1.1*font.size.scale ,
                main = "")
        )

        if(! is.null(BARCOLkinases)) {
                
                for (kin in BARCOLkinases[[i]]) {
                        
                    if(is.na(kin)) next
    
                    sel       <- names(p) %in% kin
                    q         <- unname(p)
                    q [! sel] <- 0
                    
                    idx <- which(BARCOLkinases[[i]] == kin)
    
                    barplot(q, horiz = TRUE, space = 0, yaxs = "i", axes = FALSE,
                            las = 1, col = clb[[i]][idx], main = "", add = TRUE)
                }
        }
        
        M  <- max(p)
        
        axis(side = 1, at = c(0, M), labels = FALSE, tcl = 0.2)
        
        ptck <- pretty(c(0, M), n = 4); ptck <- ptck [ptck <= M]
        
        axis(side = 1, at = ptck, padj = -0.5, tcl = -0.3, cex.axis = 1.2*font.size.scale)

        mtext("Kinase-Centric Counts", side = 1, at = 0.5*M, line = 2, cex = 1.3*font.size.scale)
        
        mtext("Not in PSP/NWK", side = 3, line = 0.6, cex = 1.3*font.size.scale, font = 3)
            
#_______________________________________________________________________________
# INKA SCATTER PLOTTING

        par(fig=c(frame.div,frame[2],frame[3],frame[4]),new=TRUE)
        par(mar=c(2,2,1,2))

        col.kin <- COLkinases[[i]]
   
        pdt [, `:=` (Kin = Kinome + ActLoop, Sub = PSP + NWK) ]
        
        pdt [, Score      := sqrt(Kin * Sub)] 
        pdt [, Rel.Score  := Score / max(Score)] 
        pdt [, Skew       := atan2(Sub, Kin) / (0.5*pi)]

        mid.position <- 0.5
        if (tiox.flag){
            cat("Centering skews in inka plot.\n")
            sel.non.zero <- pdt$Score > 0
            ratios <-  pdt$Sub[sel.non.zero] / pdt$Kin[sel.non.zero]
            med.ratio <- median(ratios)
            pdt$Skew[sel.non.zero] <- atan(ratios/med.ratio)  / (0.5*pi)
            mid.position <- atan( 1/med.ratio ) / (0.5*pi)

        }
        
        pdt [, `:=` (Multi = .N, Occur = 1:.N) , by = paste(Score, Skew)]

        pdt[, plot(Skew, Rel.Score, 
                   
                   xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n",
                   xlab = "", ylab = "", 
                   
                   pch      = 19, 
                   cex      = 1.5*font.size.scale,                                                                   
                   col      = ifelse(Kinase %in% col.kin, "red", "blue"), 
                   xlim     = c(-0.1, 1.1),                                  
                   ylim     = c(0.1,1.05), 
                   las      = 1,
                   main     = paste("INKA", sample.name),                                  
                   cex.main = 1.7*font.size.scale)                                           
        ]
        
        abline(v=c(0, 1),  lty = c("13","13"), col = "grey50")                 

        if (tiox.flag){
            axis(side   = 1, at=c(0, 1.0), labels=c("Kinase-Centric\nEvidence Only", 
                        "Substrate-Centric\nEvidence Only"),
                         padj=0.5, cex.axis=1.3*font.size.scale)
            text(mid.position,0.5,"Equal contribution",col="grey50",cex=0.75,srt=90)
            segments(mid.position,0,mid.position,0.3,col="grey50",lty="11")
            segments(mid.position,0.7,mid.position,1.1,col="grey50",lty="11")
        } else {
            axis(side   = 1, at=c(0, mid.position, 1.0), labels=c("Kinase-Centric\nEvidence Only", 
                        "Equal\nContribution","Substrate-Centric\nEvidence Only"),
                         padj=0.5, cex.axis=1.3*font.size.scale)
            abline(v=0.5,  lty = "13", col = "grey50")                 
        }
        
        axis(side=2, at=seq(0.1, 1, l=10), labels = 1:10/10, las = 1, cex.axis = 1.3*font.size.scale)

        mtext("Relative Score", side = 2, line = 2.2, cex = 1.5*font.size.scale)                      
        
        pdt[, axis(side   = 4, 
                   at     = seq(0.1, 1, l=10),  
                   labels = round(max(Score) * seq(0.1, 1, l=10), 0),
                   las = 1, cex.axis = 1.3*font.size.scale)
        ]

        mtext("Score", side = 4, line = 2, cex = 1.5*font.size.scale)                     
        
        for (idx in 1:4) {
                
            if(! idx %in% pdt$Occur) next else {                      

                pdt [Occur==idx] [, 
                              
                  text(x      = Skew, 
                       y      = Rel.Score, 
                       labels = Kinase,
                       pos = c(2, 4, 3, 1) [idx], 
                       cex =  1.2*font.size.scale,                            
                       col = ifelse(Kinase %in% col.kin, "red", cl [Multi])) 
                ]
            }
        }
        
        dt <- pdt [, .SD, .SDcols = 1:10]

        cat("Finished INKA plotting.\n")

        return(dt)
        
}

make_inka_barplot <- function(in.data, pname, frame) {
        
    font.scale <- 0.75

    cat("Making InKA barplot...")
    par(fig=frame,new=TRUE)
    dt <- copy(in.data)
    dt$Kin <- dt$KINOME + dt$ACTLOOP 
    dt$Sub <- dt$PSP + dt$NWK
    dt$Score <- sqrt(dt$Kin * dt$Sub)
    setorder(dt, -Score)
    sdt  <- dt [20:1]
    MAX  <- max(sdt$Score)

    suppressWarnings(
    sdt[, barplot(Score, names.arg = KINASE, cex.names=font.scale,
                  horiz = TRUE, space = 0, las = 1, line = -0.6,
                  axes = FALSE, yaxs = "i", col = "grey95", main = "")
             ]
    )
    axis(side = 1, at = floor(c(0, MAX)), labels = FALSE, tcl = 0.25)                
    axis(side = 1, padj = -0.5, cex.axis = 1.1*font.scale) 
    mtext("INKA Score", side = 1, line = 2.25, cex = 1*font.scale)
    mtext(paste("InKA Score ",pname), side = 3, font = 2, line = 0.8, 
            cex = 1.2*font.scale)
    cat("Ready.\n")
        
}

plot_score <- function(dt, pname, ckinases, kcl) {
        
    setorder(dt, -Score)
    sdt  <- dt [20:1]
    MAX  <- max(sdt$Score)

    suppressWarnings(
    sdt[, barplot(Score, names.arg = Kinase, 
                  horiz = TRUE, space = 0, las = 1, line = -0.6,
                  axes = FALSE, yaxs = "i", col = "grey95", main = "")
             ]
    )
    axis(side = 1, at = floor(c(0, MAX)), labels = FALSE, tcl = 0.25)                
    axis(side = 1, padj = -0.5, cex.axis = 1.1) 

    mtext("INKA Score", 
          side = 1, line = 2.25, cex = 1)
    
    mtext(pname, 
          side = 3, font = 2, line = 0.8, cex = 1.2)
        
    if (! is.null(ckinases)) {
            
        for (kin in ckinases) {
                
        if(is.na(kin)) next
            
        cdt <- sdt [, .(Kinase, Score)]
        cdt [! Kinase %in% kin, Score := 0]
        idx <- which(ckinases == kin)
        cdt [, barplot(Score,
                horiz = TRUE, space = 0, las = 1, axes = FALSE, 
                yaxs = "i", col = kcl[idx], add = TRUE)
             ]
        }
    }
}


# FUNCTION generate_inka_score_plots ###########################################

inka_score_plots <- function( plot.data,
                                      COLkinases    = NULL,
                                      cl            = NULL,
                                      panel.rows    = 2,
                                      panel.columns = 2) {
        

        cat("Starting INKA score plotting...")
                
        for(i in 1:length(plot.data)) {
                                        
            plot_score(plot.data [[i]], 
                       names(plot.data) [i],
                       COLkinases [[i]], 
                       cl [[i]])
        }
        
        cat("Ready.\n")
        
#_______________________________________________________________________________
        
}
        

#
# ------- #

permute_non_zero_spectral_counts <- function ( orig.path,out.path ){

    pp <- read.csv2(orig.path,stringsAsFactors=FALSE)
    quote.cols <- c("Proteins","Gene.Names.Original.","Gene.Names","Protein.Names","Phospho..STY..site.IDs")
    colnames.pp <- colnames(pp)
    quote.col.nrs <- match(quote.cols,colnames.pp)
    
    sample.cols <- grep("Spectral.Count.",colnames.pp)
    
    for (i in sample.cols){
        original.vals <- pp[,i]
        all.nonzero <- original.vals > 0
        values.selected <- original.vals[all.nonzero]
        selected.length <- length(values.selected)
        perm.used <- sample(selected.length)
        values.permuted <- values.selected[perm.used]
        original.vals[all.nonzero] <- values.permuted
        cat("Randomizing ", selected.length," non-zero counts of column ",colnames.pp[i],".\n")
        pp[,i] <- original.vals
    }
    write.table(pp,file=out.path,sep=";",row.names=FALSE,quote=quote.col.nrs)

}

#
# -------
#

map_gene_to_previous <- function(genes,hgnc.full){

    mappable.flag <- genes %in% hgnc.full$prev_symbol
    if( sum(mappable.flag)>0){
        genes[mappable.flag] <- hgnc.full$prev_symbol[match( genes[mappable.flag] ,hgnc.full$prev_symbol)]
    }
    return(genes)

}

#
# -------
#

gene2hgnc <- function ( genes.in ,hgnc.full,unmapped.gene="Unmapped_Gene",sloppy=FALSE){
    
    mapped.genes <- rep(unmapped.gene,length(genes.in))
    uppercased.genes <- toupper(genes.in)
    hgnc.genes <- unique(hgnc.full$symbol[ hgnc.full$status == "Approved" ])
    approved.flag <- uppercased.genes %in% hgnc.genes
    mapped.genes[approved.flag] <- uppercased.genes[approved.flag]
    unmapped.flag <- mapped.genes == unmapped.gene
    if (sum(unmapped.flag) > 0){
        cat("Unable to map right-away",sum(unmapped.flag)," genes.\n")
        mapped.genes[unmapped.flag] <- map_gene_to_previous(mapped.genes[unmapped.flag],hgnc.full)
        unmapped.flag <- mapped.genes == unmapped.gene
        if (sum(unmapped.flag) > 0){
            cat("Unable to map at all ",sum(unmapped.flag)," genes.\n")
            write.table(unique(genes.in[unmapped.flag]),file="unmapped_genes.txt",quote=FALSE,row.names=FALSE,
                    col.names=FALSE)
            if(sloppy){
                mapped.genes [unmapped.flag] <- paste0(genes.in[unmapped.flag],"*")
            }
        }
    }

    return(mapped.genes)
}

#
# -------
#

convert_to_valid_geneWise <- function ( gw.orig ,hgnc.full){

    default.gene   <- "Unmapped_Gene"

    tomap.names    <- toupper(gw.orig$Gene)
    n.gene.entries <- length(tomap.names)
    uniq.tomap     <- unique(tomap.names)
    cat("Mapping ",length(uniq.tomap)," genes: ")
    
    hgnc.approved.genes <- unique(hgnc.full$symbol[ hgnc.full$status == "Approved" ])
    genes.approved      <- uniq.tomap %in% hgnc.approved.genes
    cat(sum(genes.approved)," in HGNC, ")
    
    not.in.approved <- unique(uniq.tomap[! genes.approved])
    
    prev.symbols <- not.in.approved[not.in.approved %in% hgnc.full$prev_symbol]
    
    cat(length(prev.symbols) ," previous, ",
            sum(! not.in.approved %in% hgnc.full$prev_symbol)," unmapped. ")
    # Jaco adds these columns
    # =======================
    #
    # "Gene.Symbol.HGNC21Apr2016"    <--     FIXME only doing this one
    # "Protein.Name.HGNC21Apr2016"                     
    # "Chrom.Loc.HGNC21Apr2016"
    #
    gw.orig$Gene.Symbol.HGNC21Apr2016 <- rep(default.gene,n.gene.entries)
    
    approved.idxs <- which(tomap.names %in% hgnc.approved.genes)
    gw.orig$Gene.Symbol.HGNC21Apr2016[approved.idxs] <- gw.orig$Gene[approved.idxs]
    
    rescued <- 0
    for (i in prev.symbols){
        selected <- hgnc.full$prev_symbol == i
        hgnc.symbol <- hgnc.full$symbol[selected]
        if (hgnc.full$status[selected] == "Approved" ){
            target <- grep(i,gw.orig$Gene)
            gw.orig$Gene.Symbol.HGNC21Apr2016[target] <- hgnc.symbol
            rescued <- rescued + 1
        }
    }
    cat("[",rescued," rescued].  Ready.\n")
    return(gw.orig)
}

#
# ------------------------------------------------------------------------------------
#

add_results <- function ( big,input.file,score.name ){

    cat("Appending results from file: ",basename(input.file),"\n")
    results <- read.delim(input.file,header=TRUE,sep='\t',stringsAsFactors=FALSE)
    input.frame <- data.frame(results$Kinase,results$Score)
    colnames(input.frame) <- c("Kinase",score.name)
    if(!is.null(ncol(big)))
        big <- merge(big,input.frame,by="Kinase",all=TRUE)
    else
        big <- input.frame

    return(big)
}

#
# ------------------------------------------------------------------------------------
#

substrate_preprocess <- function ( tab ){

    tab <- tab[! duplicated(paste(Ppeptide.Id, Position.In.Peptide, KINASE))]
    tab[, N.Sites := .N, by = .(Ppeptide.Id, KINASE, SUBSTRATE)]
}

#
# ------------------------------------------------------------------------------------
#

extract_quant_col_names <- function(kin,dtype){
    candidates <- colnames(kin)[grep(dtype,colnames(kin))]
}

#
# ------------------------------------------------------------------------------------
#

construct_inka_table <- function (kin,act,psp,nwk,qname) {

    kin$signal.kin <- kin[[qname]]*kin$Phospho.Mods
    kin <- kin[,c("KINASE","signal.kin")]
    act$signal.act <- act[[qname]]*act$Phospho.Mods
    act <- act[,c("KINASE","signal.act")]
    psp$signal.psp <- psp[[qname]]*psp$Phospho.Mods/psp$N.Sites
    psp <- psp[,c("KINASE","signal.psp")]
    nwk$signal.nwk <- nwk[[qname]]*nwk$Phospho.Mods/nwk$N.Sites
    nwk <- nwk[,c("KINASE","signal.nwk")]
    cat("here we start merging stuff..\n")
    inka.kin <- merge(kin,act,by="KINASE")
    inka.sub <- merge(nwk,psp,by="KINASE")
    inka <- merge(inka.kin,inka.sub,by="KINASE")
    inka <- inka[,c("KINASE","signal.kin","signal.act","signal.psp","signal.nwk")]
    colnames(inka) <- c("KINASE","KIN","ACT","PSP","NWK")
    return (inka)
}

#
# ------------------------------------------------------------------------------------
#

extract_inka_samples <- function (kin,act,psp,nwk,dtype){

    psp <- substrate_preprocess(psp)
    nwk <- substrate_preprocess(nwk)
    sample.columns <- extract_quant_col_names(kin,dtype)
    inka.list <- list()
    for (i in seq_along(sample.columns)){
        inka <- construct_inka_table(kin,act,psp,nwk,sample.columns[i])
        inka.list[[i]] <- inka
    }
    names(inka.list) <- sub(dtype,"",sample.columns)
}

#
# ------------------------------------------------------------------------------------
#

extract_inka_scores <- function( plot.data, plot.names, dtype ) {
        
        cat("\nStarting INKA score extraction ...\n")
        
        plot.data <- lapply(plot.data, function(x) x [, sum(SIGNAL), by=KINASE])

        if(is.null(plot.names)) {
                
            cn <- nchar(dtype) + 2
            nn <- sapply(names(plot.data), 
                         function(x) substr(x, cn, regexpr("_", x ) - 1))
            
            plot.names <- unique(nn)

        }
        
        returns <- list()

        for (i in seq_along(plot.names)) {
                
            pdata   <- plot.data [1:4 + (i-1)*4]
            kins    <- unique(unlist(lapply(pdata, function(x) x [, KINASE])))
            pdt     <- data.table(KINASE = kins)
            pdt <- Reduce(function(x,y) merge(x, y, by = "KINASE", all = TRUE), pdata)
            pdt [is.na(pdt)] <- 0               
            setnames(pdt, c("Kinase", "Kinome", "ActLoop", "PSP", "NWK"))
            pdt [, `:=` (Kin = Kinome + ActLoop, Sub = PSP + NWK) ]
            pdt [, Score      := sqrt(Kin * Sub)] 
            pdt [, Rel.Score  := Score / max(Score)] 
            pdt [, Skew       := atan2(Sub, Kin) / (0.5*pi)]
            dt <- pdt [, .SD, .SDcols = 1:10]
            pn <- plot.names [i]
            returns <- c(returns, list(dt))
            names(returns)[i] <- pn
        }
        
        cat("Finished INKA score extraction.\n")
        
        return(returns)
}

#
# ------------------------------------------------------------------------------------
#

append_inka_results <- function ( seq , inka ,add.name){

    dt.tmp <- data.frame(inka$Kinase,inka$Score)
    colnames(dt.tmp) <- c("Kinase",add.name)
    merge(seq,dt.tmp,by="Kinase",all=TRUE)
}

#
# ------------------------------------------------------------------------------------
#

compose_pp_hat <- function (pp.all) {

    rep(pp.all$ppModPeptide.ID,pp.all$Global.Counts)

}

#
# ------------------------------------------------------------------------------------
#
#

fabricate_inka_inputs <- function (pp.orig,ps.orig,pp.hat,pp.all,ps.all,
                                   nonzero.spectral.counts,dtype) {

#   We select a bunch of peptides
    pp.nrows <- nrow(pp.orig)
    pp.idxs <- sample(pp.hat,pp.nrows)
    pp.columns <- c("ppModPeptide.ID","Sequence","Proteins","Gene.names",
                    "Phospho..STY.","Phospho..STY..site.IDs")
    pp <- pp.all[match(pp.idxs,pp.all$ppModPeptide.ID),pp.columns]
#   Now we add replace non-zero sample counts
    sample.col.idxs <- grep(dtype,colnames(pp.orig))
    sample.col.names <- colnames(pp.orig)[sample.col.idxs]
    n.samples <- length(sample.col.idxs)
    for (i in seq_along(sample.col.idxs)){
        nz.idxs <- pp.orig[[ sample.col.idxs[i] ]] > 0
        new.vals <- sample(nonzero.spectral.counts, sum(nz.idxs) )
        new.df <- data.frame(counts=rep(0,pp.nrows))
        new.df$counts [nz.idxs] <- new.vals
        colnames(new.df) <- sample.col.names[i]
        pp <- cbind(pp,new.df)
    }
    tmp.names <- colnames(pp)
    tmp.names <- sub("Gene.names","Gene.Names",tmp.names)
    colnames(pp) <- tmp.names

#   Construct the corresponding Phospho (STY)Sites.txt file
    ps.ids <- unique(unlist(strsplit(pp$Phospho..STY..site.IDs,";")))
    ps.colnames <- c("Proteins","Positions within proteins","Localization prob",
            "Amino acid","Sequence window","Reverse","Contaminant","id","Mod. peptide IDs")
    ps <- ps.all[match(ps.ids,ps.all$id),]
    colnames(ps) <- ps.colnames
#   This is what we return
    retval <- list(data.table(pp),data.table(ps))
    names(retval) <- c("PP","PS")
    return (retval)

}

#
# ------------------------------------------------------------------------------------
#

inka_plotdata_mangle <- function (kin,act,psp,nwk,dtype,sample.names) {

    foldSort( relabel( c(getData(kin,dtype), getData(act,dtype), getData2(psp,dtype), getData2(nwk,dtype)), 
                       sample.names,dtype) )

}

#
# ------------------------------------------------------------------------------------
#

find_sample_data_indices <- function (raw,dtype,sample.name) {

    object.names <- names(raw)
    kin.tag <- paste0(dtype,".",sample.name,"_Kinome")
    act.tag <- paste0(dtype,".",sample.name,"_ActivationLoop")
    psp.tag <- paste0(dtype,".",sample.name,"_PhosphoSitePlus")
    nwk.tag <- paste0(dtype,".",sample.name,"_NetworKIN")
    kin.idx <- match(kin.tag,object.names)
    act.idx <- match(act.tag,object.names)
    psp.idx <- match(psp.tag,object.names)
    nwk.idx <- match(nwk.tag,object.names)
    retval <- c(kin.idx,act.idx,psp.idx,nwk.idx)

    return(retval)


}

#
# ------------------------------------------------------------------------------------------------------
#

extract_sample_dataset <- function ( raw,dtype,sname ){

    data.idxs <- find_sample_data_indices(raw,dtype,sname)
    retval <- raw[data.idxs]

    return(retval)
}

#
# ------------------------------------------------------------------------------------------------------
#

opl_pp_read <- function ( filename,dtype=NULL,labfile="labels.txt" ){

    if(is.null(dtype)){
        dtype <- "Spectral.Counts"
    }
    raw <- fread(filename, integer64 = "numeric")
    raw$Gene.names <- gsub("\"","",sub("=","",raw$Gene.names))
    colnames(raw) <- sub("Gene.names","Gene.Names",colnames(raw))
    samples <- grep(dtype,colnames(raw))
    if (length(samples) == 0){
        colnames(raw) <- sub("OPL.ppCount.",paste0(dtype,"."),colnames(raw))
    }

    if(file.exists(labfile)){
        cat("Found label file. Fixing names.\n")
        labs <- read.delim(labfile,sep="\t",stringsAsFactors=FALSE,header=FALSE)
        # Get rid of obnoxious names
        labs[[ncol(labs)]] <- sub("\\.+$","",sub("^\\.+","",make.names(labs[[ncol(labs)]])))
        colnames(raw) <- opl_fix_colnames(raw,dtype,labs)
    }

    return(raw)
}

#
# ------------------------------------------------------------------------------------------------------
#

opl_fix_colnames <-function ( pp, dtype, lab.defs ){

    old.names <- colnames(pp)
    last.col <- ncol(lab.defs)
    new.names <- make.names(lab.defs[,last.col])
    for (i in 1:nrow(lab.defs)){
        in.string  <- paste0(dtype,".",lab.defs[i,1],"$")
        out.string <- paste0(dtype,".",new.names[i])
        old.names  <- sub(in.string,out.string,old.names)
    }
    return(old.names)
}

#
# ------------------------------------------------------------------------------------------------------
#



trimComponents <- function(edges) { 
               
    cln1_2 <- names(edges) [1:2] 
    EE     <- copy(edges)
    setnames(EE, c("V1", "V2", names(EE)[-(1:2)]))
            
    vertices <- EE [, .(VERTEX = unique(c(V1, V2)))]
    
    VV <- copy(vertices)
    VV [, `:=` (Discovered = FALSE, Processed = FALSE)]

    component <- 0
    proceed   <- TRUE
 
    while(proceed) {                  

        component <- component + 1
    
        VV [1, Discovered := TRUE]
            
        tagger <- VV [1, VERTEX]
        tagged <- EE [grep(tagger, Vertex.Pair), unique(c(V1,V2))] 
        VV [VERTEX %in% tagged, Discovered := TRUE]
            
        VV [1, Processed := TRUE]
        
        continue <- TRUE
        
        while(continue) {
                
            taggers <- VV [Discovered & !Processed, VERTEX]
                 
            if(length(taggers) > 0) {
                 
                for(tagger in taggers){
                         
                    tagged <- EE [grep(tagger, Vertex.Pair), unique(c(V1,V2))]
                    VV [VERTEX %in% tagged, Discovered := TRUE]
                    VV [VERTEX == tagger, Processed := TRUE]
                }
                                 
            } else continue <- FALSE
        }
    
        marked <- VV [{Discovered}, VERTEX]
        vertices [VERTEX %in% marked, sprintf("Component%s", component) := 1]
        
        VV <- VV [{! Discovered}]
        
        if(nrow(VV) == 0) proceed <- FALSE
    }
    
    vertices[is.na(vertices)] <- 0
    tmp <- names(vertices)
    setcolorder(vertices, c(1, order(colSums(vertices[, -1]), decreasing = TRUE) + 1))
    setnames(vertices, tmp)
    
    major.vertices  <- vertices[Component1 == 1, VERTEX]
    EE <- EE [V1 %in% major.vertices]
    setnames(EE, c(cln1_2, names(EE)[-(1:2)]))

    #commented out: return(list(trimmed.Edges = EE, Edges = edges, Vertices = vertices))
    return(EE)
}

# ----------------------------

prepare_KSR <-function (psp_ps_tab,nwk_ps_tab,actloop,kinase.symbols,dtype,hgnc.full ){

    colnames(psp_ps_tab) <- make.names(colnames(psp_ps_tab))
    colnames(nwk_ps_tab) <- make.names(colnames(nwk_ps_tab))

    KSR <- merge(x   = psp_ps_tab,
                 y   = nwk_ps_tab,
                 by  = c("KINASE", "SUBSTRATE", "Ppeptide.Id", "Psite.Id", "Accession_AApos"), all = TRUE)
    
    KSR [, SUBSTRATE := gene2hgnc(SUBSTRATE,hgnc.full,sloppy=TRUE)]
    
    dtnames <- names(psp_ps_tab)
    snames  <- dtnames [grep(dtype, dtnames)]    
    snames  <- sapply(snames, function(x) str_sub(x, start = nchar(dtype)+2))

    KSR [, `:=` ( Source      = mapply(function(x,y) { if(is.na(x)) "NWK" else if(is.na(y)) "PSP" else "Both" },
                                      Sequence.x, Sequence.y),  
                  Vertex.Pair = mapply(function(x,y) { paste(sort(c(x, y)), collapse = " ") }, KINASE, SUBSTRATE),
                  Sub.is.Kin  = SUBSTRATE %in% kinase.symbols,
                  Auto.Phos   = KINASE == SUBSTRATE,
                  Loop.Site   = Ppeptide.Id %in% actloop$Ppeptide.Id ) ]
     
    for (i in seq_along(snames)) {

        tmp <- paste0(dtype,".",snames[i])
#tmp <- snames[i]
        cat("Working on",tmp,"\n")
        clm <- grep(tmp, names(KSR))
        cat("Column names processed:",names(KSR)[clm],"\n")
        
        KSR [, paste0(tmp, ".pepSIGNAL") := rowSums(.SD, na.rm = TRUE), .SDcols = clm] 
        KSR [Source == "Both"] [, paste0(tmp,".pepSIGNAL") := eval(as.name(paste0(tmp, ".pepSIGNAL"))) / 2]
        KSR [, paste0(tmp, ".SIGNAL") := sum(eval(as.name(paste0(tmp, ".pepSIGNAL")))), by = .(KINASE, Accession_AApos)]
    }

    rmv <- unlist(sapply(c("[.](x|y)$", "P(site|peptide)[.]Id", "[.]pepSIGNAL"), function(x) grep(x, names(KSR))))
    KSR <- KSR [, .SD, .SDcols = -rmv]
    KSR <- KSR [! duplicated(paste(KINASE, Accession_AApos))]

    return(list(snames=snames,KSR=KSR))
}
# ----------------------------

select_sample_KSR <- function (KSR,snames,sname,dtype){

#cat("From all sample names: ",snames,"\n")
#cat("selecting: ",sname,"\n")
    sname.simplified <- gsub("\\.\\.*","\\.",sname)
    leaveout.names <- paste0(dtype,".",snames[snames != sname.simplified],".SIGNAL")
#cat("Leave-out name: ",leaveout.names,"\n")
    clm   <- match(leaveout.names, names(KSR) )
#cat("In function select_sample_KSR ignoring column :",clm,"\n")
#cat("length(clm)= ",length(clm),"\n")
#EDGES <- if(length(clm) > 0) KSR [, .SD, .SDcols = -clm] else copy(KSR)
    EDGES <- if(sum(!is.na(clm)) > 0) KSR [, .SD, .SDcols = -clm] else copy(KSR)
    col.idx <- grep ("SIGNAL",colnames(EDGES))
    colnames(EDGES)[col.idx] <- "SIGNAL"
    EDGES <- EDGES [SIGNAL > 0]

    return(EDGES)
}

#
# ----------------------------
#

generate_KSR_network <- function (EDGES,
                                  inka,
                                  sname,
                                  kinase.symbols,
                                  signal_cutoff    = NULL,
                                  topN             = NULL,
                                  major_component  = FALSE,
                                  plot_large       = TRUE,
                                  non_Kinase_sides = 50,
                                  obs_Kinase_sides = 6,
                                  inf_Kinase_sides = 4,
                                  palette_steps    = 100,
                                  pal_color_from   = "white", 
                                  pal_color_to     = "red",
                                  non_Kinase_color = "#DBDBDB",
                                  edge_colors      = c(PSP  = "coral", NWK  = "cornflowerblue", Both = "forestgreen")
                                  ) {

    cat("Creating KSR network....")

    kinase.colors <- colorRampPalette(c(pal_color_from, pal_color_to))(palette_steps + 1)
            
    setnames(inka, c("VERTEX", names(inka)[-1]))
    setorder(inka, -Score, na.last = TRUE)
    
    if(! is.null(signal_cutoff)) {
            
        cutoff.operator <- str_replace_all(signal_cutoff, "[0-9]", "")
        cutoff.value    <- as.integer(str_replace_all(signal_cutoff, "[<>=]", ""))
        operator.name   <- c("Greater","GreaterEq", "LessEq","Less")[c(">",">=", "<=", "<") %in% cutoff.operator]
        EDGES <- EDGES [mapply(cutoff.operator, SIGNAL, cutoff.value)]
    }
    
    if(! is.null(topN)) {

        topN.kinases <- inka [1:topN, VERTEX]
        EDGES <- EDGES [KINASE %in% topN.kinases | SUBSTRATE %in% topN.kinases]
    }
           
    if(major_component) {
            
        EDGES <- trimComponents(EDGES)
    }

    n.vertices <- nrow(EDGES)
    cat("Number of edges : ",n.vertices,"\n")
    if (n.vertices){
       
    EDGES [, N.Edges := .N, by = Vertex.Pair]
    EDGES [, N.Forward.Edges := .N, by = .(KINASE, SUBSTRATE)]
    EDGES [, N.Reverse.Edges := N.Edges - N.Forward.Edges]
    
    EDGES [, Curve := 0.1*(1:.N), by = .(KINASE, SUBSTRATE)]
    EDGES [, Color.Idx := match(Source, c("PSP", "NWK", "Both"))]
    EDGES [, Rel.Signal := SIGNAL / max(SIGNAL)]

    VERTICES <- EDGES [, .(VERTEX = unique(c(KINASE, SUBSTRATE)))]
    VERTICES <- merge(x  = VERTICES,
                      y  = inka,
                      by = "VERTEX", all.x = TRUE)
    
    VERTICES [, is.Kinase := VERTEX %in% kinase.symbols]
    VERTICES [, Sides := ifelse(! is.Kinase, non_Kinase_sides,
                                ifelse(Kinome > 0, obs_Kinase_sides, inf_Kinase_sides))]
    VERTICES [, Color := ifelse(! is.Kinase, non_Kinase_color,
                                kinase.colors [round(100 * Rel.Score) + 1])]
    VERTICES [, Loop.Site := VERTEX %in% EDGES[{Loop.Site}, SUBSTRATE]]

        

        VERTICES[is.na(VERTICES)] <- 0
        VERTICES$VERTEX <- as.character(VERTICES$VERTEX)
        NW <- network(EDGES,
                      vertex.attr = VERTICES,
                      loops       = TRUE,
                      multiple    = TRUE,
                      matrix.type = "edgelist",
                      ignore.eval = FALSE)
        size       <- network.size(NW)
        layout.par <- list(niter              = size * 100,
                           area               = size^1.8,
                           repulse.rad        = size^1.5,
                           ncell              = size^3)


        set.seed(1)

        plot.network(NW,
                     new           = TRUE,
                     arrowhead.cex = 5,
                     boxed.labels  = if(plot_large) FALSE else TRUE,
                     coord         = network.layout.fruchtermanreingold(NW, layout.par),
                     displaylabels = TRUE,
                     edge.col      = edge_colors [NW %e% "Color.Idx"],
                     edge.curve    = 0.75 * NW %e% "Curve",
                     edge.lty      = 1,
                     edge.lwd      = 1 + 25 * (NW %e% "Rel.Signal"),
                     edge.steps    = 50,
                     label.bg      = "#000000",
                     label.border  = "#000000",
                     label.cex     = if(plot_large) 0.5 else 0.4,
                     label.col     = if(plot_large) "#000000" else "#FFFFFF",
                     label.pad     = 0.2,
                     label.pos     = 5,
                     loop.cex      = 1,
                     loop.steps    = 20,
                     mode          = "fruchtermanreingold",
                     object.scale  = 0.001,
                     suppress.axes = TRUE,
                     usearrows     = TRUE,
                     usecurve      = TRUE,
                     vertex.border = ifelse(VERTICES$Loop.Site, "#4D4D4D", "#333333"),
                     vertex.cex    = if(plot_large) 25 else 15,
                     vertex.col    = NW %v% "Color",
                     vertex.lty    = 1,
                     vertex.lwd    = if(plot_large) ifelse(VERTICES$Loop.Site, 3, 1) else 
                                                    ifelse(VERTICES$Loop.Site, 3, 1),
                     vertex.sides  = NW %v% "Sides",
                     vertices.last = TRUE)

         mtext(sname, side = 3, cex = 1.5, line = 1, font =2)
    } else {
        textplot("Not enough data for network.")
    }
        
    cat("Ready.\n")

}


#
# ------------------------------------------------------------------------------------------------------
#

extract_mq_version <- function ( pp.file ){

    first.line <- readLines(pp.file, n=1) 
    c.names <- unlist(strsplit(first.line,split="\t"))
    found <- grep("ontaminant",c.names)
    if (length(found)==0){
        version <- "1.5"
    }
    else {
        if (c.names[found] == "Potential contaminant"){
            version <- "1.5"
        }
        else {
            version <- "1.4"
        }
    }
    return(version)
}

#
# ------------------------------------------------------------------------------------------------------
#
cat("This is inka.R version ",inka.script.version,"\n")
cat("Loading necessary data...")
UniProt_data  <- loadf(paste0(general.resources.dir,"/","UniProt_data_08Jun2016.Rdata"))
PSP_human     <- loadf(paste0(general.resources.dir,"/","PSP_human_03Jul2016.Rdata"))
PSP_KSR_human <- loadf(paste0(general.resources.dir,"/","PSP_KSR_human_03Jul2016.Rdata"))
Kinases       <- loadf(paste0(general.resources.dir,"/","Manning.Kinases_20Jan2017.Rdata"))
hgnc.resources.file <- paste0(general.resources.dir,"/","complete-HGNC-21Apr2016.csv")
hgnc.full           <- read.csv2(hgnc.resources.file,stringsAsFactors=FALSE)
phomics.path        <- paste0(general.resources.dir,"/","phomics_global_output.txt")
phomics.global      <- fread(phomics.path)
nwk.global.path     <- paste0(general.resources.dir,"/","TOP_NWK_nwk_proteome_unprot_ref_2014.Rdata")
load(nwk.global.path)
cat(" Done.\n")

mapgen     <- "Gene.Symbol.HGNC21Apr2016"
dtype      <- "Spectral.Count"

#mq.version <- "1.5"
mq.version <- extract_mq_version(ps.filepath)
cat("Detected MAxQuant version ",mq.version,"\n")

pp.orig.1CL <- opl_pp_read(pp.filepath,dtype)

samples <- sub(paste0(dtype,"."),"",colnames(pp.orig.1CL)[grep(dtype,colnames(pp.orig.1CL))])
cat("Extracted sample names:\n")
for ( i in seq_along(samples)){
    cat(" ",samples[i])
    if(i %% 10 == 0) cat ("\n")
}
cat("\n")

cat("Reading phospho-site file:",ps.filepath,"..")
ps.orig.1CL <- read_mq_ps_file(ps.filepath,mq.version)
cat(" Ready.\n")
# ---
# Filter annotation step
# output of filter_annotate: list(report_PP,geneWise_report_PP,nonRed_classI_PS_PP)


fa_results.1CL <- filter_annotate( pp.orig.1CL,ps.orig.1CL,UniProt_data,PSP_human,PSP_KSR_human,mq.version)
kin.1CL        <- generate_kinome_pdata(convert_to_valid_geneWise(fa_results.1CL[[2]],hgnc.full),Kinases,mapgen)
actloop.1CL    <- generate_act_loop_pdata(phomics.global,kin.1CL,Kinases,mapgen,dtype)
psp.1CL        <- generate_PSP_pdata(fa_results.1CL[[3]],PSP_KSR_human,mapgen,dtype)
#save(psp.1CL,file="PSP_dump.Rdata")
nwk.1CL        <- generate_NWK_pdata(TOP_NWK,fa_results.1CL[[3]],dtype,mapgen)

full.data <- inka_plotdata_mangle(kin.1CL,actloop.1CL,psp.1CL,nwk.1CL,dtype,samples)

cat("Preparing for drawing KS networks.\n")
netKSR <- prepare_KSR (psp.1CL,nwk.1CL,actloop.1CL,Kinases$Gene.Symbol.HGNC21Apr2016,dtype,hgnc.full )

cat("Writing output to PDF file ",inka.pdf.output.filename,"\n")
document.title <- paste0("InKAscore.org phosphoproteomics data analysis version ",inkascore.version,
        " running algorithm version ",inka.script.version)
pdf(file=inka.pdf.output.filename,paper="a4",width=8,height=11,
        title=document.title)

for (i in seq_along(samples)){

    par(new=FALSE)
    cat(i,") Sample ",samples[i]," being added to report.\n")
    full.sample <- extract_sample_dataset(full.data,dtype,samples[i])
    frame <- c(0.05,0.95,0.75,1.0)
    make_kinase_bar_plot(full.sample,samples[i],frame)
    inka.data <- kinase_condense_data(full.sample)
    par(new=TRUE)
    frame <- c(0.05,0.95,0.4,0.7)
    inka.save <- make_inka_plot(inka.data,samples[i],mapgen,Kinases,frame,tiox.flag=tiox.mode)
    write.table(inka.save,file=paste0("inka_save_",samples[i],".txt"),row.names=FALSE,quote=FALSE,sep="\t")
    frame <- c(0.05,0.45,0.0,0.35)
    make_inka_barplot(inka.data,samples[i],frame)
    frame <- c(0,1.0,0.0,0.9)
    par(fig=frame,new=FALSE)
    generate_KSR_network(select_sample_KSR(netKSR[["KSR"]],netKSR[["snames"]],samples[i],dtype) ,
            inka.save,samples[i],Kinases$Gene.Symbol.HGNC21Apr2016,signal_cutoff = network.link.threshold,
            topN = 20)

}

dev.off()
