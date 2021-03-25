#' Converts annotated DMRs in a plot
#'
#' This function allows to track in a plot the beta value of a methylated segment mapped on a transcript of the human genome
#' 
#' @param annotatedDMRs anotated DMRs list resultingfrom annotateDMRs() or scoreAnnotatedDMRs()
#' @param symbol character, gene symbol to plot.
#' @param annotation character, "ensembl" or "ucsc". Annotation used to track the plot. Default = "esembl".
#' @param hg character, "hg19" or "hg38". Genome Assembly version. Default = "hg19".
#' @param beta1.name character, if unsued beta difference is plotted. character string identifying beta value of first sample in "other" column in results from annotateDMRs() or indentifying colname in input table used in annotateDMRs()
#' @param beta2.name character, it identifies beta value of second sample in "other" column in results from annotateDMRs() or it identifies colname in input table used in annotateDMRs()
#' @param beta.colors character vectors, colors of tracks for the first and the second bvalue, respectivetely. Default is c("red", "navy"). If beta diff is plotted, only the first element of vector is considered.
#' @param blackandwhite logical, it allows to get all the plot in greyscale. Default = FALSE.
#' @param show.all.transcripts logical, if TRUE all transcripts of genes are tracked, if FALSE only the longest transcript is tracked. Default = FALSE.
#' @param prom.width integer, promoter lenght. Default = 1500.
#' @param path logical, path where the plot is saved in a pdf file. If NULL the plot is not saved. Default = NULL.
#' @param coord.zoom numeric vectors, coordinates of zoom region. If NULL the plot is not zoomed. Default = NULL.
#' @param smartzoom logical, automatic zoom on the methylated region. Default = TRUE.
#' @param height.pdf integer, hight pdf file. Default = 9.
#' @param width.pdf integer, width pdf file. Default = 16.
#'
#' @return Plot of the beta value(s) mapped on transcript(s).
#'
#' @export

plotDMRs = function(annotatedDMRs, symbol, annotation = "ensembl", hg = "hg19", beta1.name=NULL, beta2.name = NULL, beta.colors = c("red","navy"), blackandwhite = FALSE, show.all.transcripts = FALSE, prom.width = 1500, path = NULL, coord.zoom = NULL, smartzoom = TRUE, height.pdf = 9, width.pdf= 16){
    symbol <- toupper(symbol)
    libpath <- paste0(.libPaths()[1], "/methyl.O")
    
    if(annotation == "ucsc"){
        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
        cgis = readRDS(paste0(libpath,"/data/CGIs.19.gr.RDS")) # load all CGIs genomic granges
        tab.all.trx = readRDS(paste0(libpath, "/data/txsSelectedAll_UCSC_hg19.RDS")) # load all trx ucsc
        tab.trx.max = readRDS(paste0(libpath,"/data/txsSelectedLongest_UCSC_hg19.RDS")) #load all trx max ucsc
        if(prom.width == 1500){
            prom = readRDS(paste0(libpath,"/data/promSelectedAll_UCSC_hg19.RDS")) #load all prom 1500 upstream ucsc
        } else {
            prom = as.data.frame(GenomicRanges::promoters(txdb, upstream = prom.width, downstream = 0, columns = c("gene_id","tx_name"))) #find all prom no 1500 upstream ucsc
        }
	} else {
        if(hg == "hg19"){
            cgis = readRDS(paste0(libpath,"/data/CGIs.19.gr.RDS")) # load all CGIs genomic granges
            edb = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
            tab.allTrack.for.Gviz.Ensdb= as.data.frame(readRDS(paste0(libpath,"/data/AllTrack_for_Gviz_Ensdb75.RDS"))) #get all trx track for Gviz ensembl 75
            ensembldb::seqlevelsStyle(edb) <- "UCSC"
            tab.all.trx = readRDS(paste0(libpath, "/data/txsSelectedAll_ENSEMBL_hg19.RDS")) # load all trx ensembl 75
            tab.trx.max = readRDS(paste0(libpath, "/data/txsSelectedLongest_ENSEMBL_hg19.RDS")) # load all max trx ensembl 75
            if(prom.width == 1500){
                prom = readRDS(paste0(libpath, "/data/promSelectedAll_ENSEMBL_hg19.RDS")) # find all prom 1500 upstream ensembl 75
            } else {
                prom = suppressWarnings(as.data.frame(GenomicRanges::promoters(edb, upstream = prom.width, downstream = 0, columns = c("gene_id","tx_name", "tx_biotype")))) #find all prom no 1500 upstream ensembl 75
            }
        } else {
            cgis = readRDS(paste0(libpath,"/data/CGIs.38.gr.RDS")) # load all CGIs genomic granges
            edb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
            ensembldb::seqlevelsStyle(edb) <- "UCSC"
            tab.allTrack.for.Gviz.Ensdb = readRDS(paste0(libpath, "/data/AllTrack_for_Gviz_Ensdb86.RDS")) #get all trx track for Gviz ensembl 86
            tab.all.trx = readRDS(paste0(libpath,"/data/txsSelectedAll_ENSEMBL_hg38.RDS")) #load all trx ensembl 86
            tab.trx.max = readRDS(paste0(libpath,"/data/txsSelectedLongest_ENSEMBL_hg38.RDS")) #load all trx max ensembl 86
            if(prom.width == 1500){
                prom = readRDS(paste0(libpath,"/data/promSelectedAll_ENSEMBL_hg38.RDS")) #find all prom 1500 upstream ensembl 86
            } else {
                prom = suppressWarnings(as.data.frame(promoters(edb, upstream = prom.width, downstream = 0, columns = c("gene_id","tx_name", "tx_biotype")))) #find all prom no 1500 upstream ensembl 86
            }
        }
	}


	#filter input annotatedDMRs[1] through symbol gene
	tab.temp <- as.data.frame(annotatedDMRs[[1]])
	tab.temp <- tab.temp[which(tab.temp$symbol %in% symbol), ]
	ind.rmv <- which(duplicated(tab.temp[, c(1:3)]))
	if (length(ind.rmv)) {
    	tab.temp <- tab.temp[-ind.rmv, ] # remove duplicated methylated segment
	}
    col.hold <- c('seqnames','start','end','width','beta','gene.start','gene.end','gene.width','gene.strand','gene.id','tx.name','others','symbol')
    tab.temp = tab.temp[, col.hold] 

	#extract the 2 bvalues in columns "others" and add this as new 2 columns to tab_temp
	col.others = as.character(tab.temp[["others"]])
	beta1 = c()
	beta2 = c()
	if (!is.null(beta1.name)) {
		for (i in 1:length(col.others)) {
			col.others.cur <- as.character(col.others[i])
			tag.beta1 <- regmatches(col.others.cur, regexpr(paste0(beta1.name, "=\\d\\.\\d+"), col.others.cur))
			bvalue1.cur <- abs(as.numeric(gsub("[^ ]+=", "", tag.beta1)))
			beta1 <- c(beta1, bvalue1.cur)
			beta1 <- as.data.frame(beta1)
			tab.temp = cbind(tab.temp[, c(1:5)], beta1, tab.temp[, c(6:13)])
		}

	if (!is.null(beta2.name)) {
		for (i in 1:length(col.others)) {
			col.others.cur <- as.character(col.others[i])
			tag.beta2 <- regmatches(col.others.cur, regexpr(paste0(beta2.name, "=\\d\\.\\d+"), col.others.cur))
			bvalue2.cur <- abs(as.numeric(gsub("[^ ]+=", "", tag.beta2)))
			beta2 <- c(beta2, bvalue2.cur)
		}
		tab.temp = cbind(tab.temp[, c(1:5)], beta1, beta2, tab.temp[, c(7:14)])
	}
	} else {
	colnames(tab.temp)[5] <- "beta1"
	}
	chr.gene =	as.character(unique(tab.temp[,"seqnames"]))
	str.met = tab.temp[,"start"] #str meth region
	end.met = tab.temp[,"end"] #end meth region

	#set the width track of meth segment: in horizontal (if only one) or vertical (more than one)
	end.start = c()
	for(i in 1:nrow(tab.temp)){
        end.start.cur = c(tab.temp[,"end"][i] - tab.temp[,"start"][i+1])
        end.start = c(end.start,end.start.cur)
	}
	end.start = abs(end.start[-which(is.na(end.start))])
	if(length(which(end.start < 5000)) >= 1){
	rotation.item = 90
	} else {
    	rotation.item = 0
	}

	#set grafic name parameters of track (name of track, on the left side of the plot)
	if(show.all.transcripts == T){
        rotation.title = 0
        title.width = 2.3
        title.width.save = 2.3
        name1 = "                                                b Value"
        if(!is.null(beta1.name) & is.null(beta2.name)){
          name2 = paste0("                                                       \"",beta1.name,"\"",sep = " ","Value")
        }else if(!is.null(beta2.name) & is.null(beta1.name)){
          name2 = paste0("                                                       \"",beta2.name,"\"",sep = " ","Value")
        }else{
          name2 = paste0("                                                       \"","Beta","\"",sep = " ","Value")
        }
        cex.axis = 0.7
        title = paste0(symbol, sep = "_", "all_transcripts")
        title.zoom = paste0(title, sep ="_", "zoom")
        name.width = "											Width"
        cex.title.metTrack = 0.9
	} else {
        rotation.title = 90
        title.width = 1.9
        title.width.save = 1.9
        name1 = "b Value"
        if(!is.null(beta1.name) & is.null(beta2.name)){
          name2 = paste0("\"",beta1.name,"\"",sep = " ","Value")
        }else if(!is.null(beta2.name) & is.null(beta1.name)){
          name2 = paste0("\"",beta2.name,"\"",sep = " ","Value")
        }else{
          name2 = paste0("\"","Beta","\"",sep = " ","Value")
        }
        cex.axis = 0.8
        title = paste0(symbol, sep = "_", "transcript_max")
        title.zoom = paste0(title, sep ="_", "zoom")
        name.width = "Width"
        cex.title.metTrack = 1
	}
	if(blackandwhite){
        background.title <- "gray86"
        background.panel <- "white"
        beta.colors <- c("gray23", "gray53")
	}else{
        background.title = "firebrick"
        background.panel = "whitesmoke"
	}


	#if there is only one bvalue to plot type = histogram, if there are two bvalues to plot type = s (line, like histogram)
	if(is.null(beta2.name)){
        type = "histogram"
        title = paste0(title,sep = "_",colnames(tab.temp$beta1))
        title.zoom = paste0(title,sep = "_",colnames(tab.temp$beta1))
	} else {
        type = "s"
        title = title
        title.zoom = title.zoom
	}

	#if there are coord.zoom, define these in function parameters: from.zoom and to.zoom
	if(!is.null(coord.zoom)){
        from.zoom = coord.zoom[1]
        to.zoom = coord.zoom[2]
	} else {
        coord.zoom = NULL
    }

	#overlap genomic region: CGIs - trx
	trxs.genes.int = tab.all.trx[which(tab.all.trx[,"gene.id"] == as.character(unique(tab.temp[,"gene.id"]))),]
	trxs.genes.int = trxs.genes.int[which(trxs.genes.int$seqnames %in% chr.gene),]
	trx.max.genes.int = trxs.genes.int[which.max(trxs.genes.int$width),]
	prom.trxs.genes.int = as.data.frame(prom[which(prom$tx_name %in% trxs.genes.int$tx.name),])
	prom.trx.max.int = prom.trxs.genes.int[which(prom.trxs.genes.int[["tx_name"]] %in%	trx.max.genes.int[["tx.name"]]),]
	strand = as.character(trx.max.genes.int[,"strand"])
	str.gene = trxs.genes.int[which.min(trxs.genes.int$start),][["start"]]
	end.gene = trxs.genes.int[which.max(trxs.genes.int$end),][["end"]]
	ind.gene.cgis = data.frame("seqnames" = chr.gene, "start" = str.gene, "end" = end.gene)
	hits = GenomicRanges::findOverlaps(cgis, GenomicRanges::makeGRangesFromDataFrame(ind.gene.cgis))
	overlapping.cgis = data.frame(cgis[S4Vectors::queryHits(hits), ])
	str.cgis.gene = overlapping.cgis$start
	end.cgis.gene = overlapping.cgis$end
	ranges.CGIs = data.frame(seqnames = rep(chr.gene,length(str.cgis.gene)), start.cgis.gene = str.cgis.gene, end.cgis.gene = end.cgis.gene)

	#get the track with overlaps CGIs:trx
	if(nrow(ranges.CGIs) == 0){
	    cgis.gene =	GenomicRanges::makeGRangesFromDataFrame(data.frame(chr = chr.gene,	start = 100000, end = 200000))
	} else {
        cgis.gene =	GenomicRanges::makeGRangesFromDataFrame(data.frame(chr = chr.gene,	start = str.cgis.gene, end = end.cgis.gene))
	}

	#get from (= start trx) and to (= end trx) for the genomictrack. FC is a correction factor releted to gene.width
	if(unique(tab.temp[,"gene.width"] <= 1000)){
	    FC = 100
	}else if(unique(tab.temp[,"gene.width"] > 1000 & unique(tab.temp[,"gene.width"] <= 5000))){
	    FC = 500
	}else if(unique(tab.temp[,"gene.width"] > 5000 & unique(tab.temp[,"gene.width"] <= 10000))){
	    FC = 1000
	}else if(unique(tab.temp[,"gene.width"] > 10000 & unique(tab.temp[,"gene.width"] <= 20000))){
	    FC = 5000
	} else {
	    FC = 10000
	}
	if(show.all.transcripts == T){
        str.prom = unique(prom.trxs.genes.int[[2]]) #info for promoter track
        if(strand == "-"){
            from = trxs.genes.int[which.min(trxs.genes.int[,"start"]),"start"] - FC
            to = prom.trxs.genes.int[which.max(prom.trxs.genes.int[,"end"]),"end"] + FC
        } else {
            from = prom.trxs.genes.int[which.min(prom.trxs.genes.int[,"start"]),"start"] - FC
            to = trxs.genes.int[which.max(trxs.genes.int[,"end"]),"end"] + FC
        }
	} else {
        str.prom =	prom.trx.max.int[["start"]] #info for promoter track
        if(strand == "-"){
            from = trx.max.genes.int[,"start"] - FC
            to = prom.trx.max.int[,"end"] + FC
        } else {
            from =	prom.trx.max.int[,"start"] - FC
            to = trx.max.genes.int[,"end"] + FC
        }
	}

	###create all track (6) of genomic info (Chr info, CGIs info, Promoter, track trx/trxs, bvalue/bvalues, width meth segment)

	#track for genomic coordinates info (axis and chromosome)
	## if black and white is TRUE
    if (blackandwhite) {
        track.gene.genome.axis <- Gviz::GenomeAxisTrack(add53 = TRUE, add35 = TRUE, showTitle = T, name = "Chr Info", rotation.title = T, fontcolor = "black", cex.title = 3, distFromAxis = 0.5, min.width = 0.5, fontsize = 5)
        track.gene.chromosome <- Gviz::IdeogramTrack(genome = hg, chromosome = chr.gene, fontsize = 5, fontcolor = "black", col = "gray73", fill = "gray73")
    } else {
        track.gene.genome.axis <- Gviz::GenomeAxisTrack(add53 = TRUE, add35 = TRUE, showTitle = T, name = "Chr Info", rotation.title = T, fontcolor = "black", cex.title = 3, distFromAxis = 0.5, min.width = 0.5, fontsize = 5)
        track.gene.chromosome <- Gviz::IdeogramTrack(genome = hg, chromosome = chr.gene, fontsize = 5, fontcolor = "black")
    }

	#track for CGIs info
	if(nrow(ranges.CGIs) == 0){
    	track.gene.CpG.island = Gviz::AnnotationTrack(cgis.gene, genome = hg, name = "CGIs site", fill= "whitesmoke", fontsize = 20, rotation.title = 0, cex.title = 0.7)
	}else{
	    if(blackandwhite){
           track.gene.CpG.island <- Gviz::AnnotationTrack(cgis.gene, genome = hg, name = "CGIs site", fill = "gray23", fontsize = 20, rotation.title = 0, cex.title = 0.7)
        }else{
           track.gene.CpG.island <- Gviz::AnnotationTrack(cgis.gene, genome = hg, name = "CGIs site", col = "purple", fill = "purple", fontsize = 20, rotation.title = 0, cex.title = 0.7)
	    }
	}

	#track for promoter
	if(strand == "-"){
	  direction.arrow = "<-"
	}else{
	  direction.arrow = "->"
	}
	direction.arrow = rep(direction.arrow, length(str.prom))
    if (blackandwhite) {
        track.promoter.max.trx <- Gviz::AnnotationTrack(start = str.prom, width = prom.width, chr = chr.gene, strand = strand, hg = hg, name = "Promoter", id = direction.arrow, showTitle = T, fill = "gray73", shape = "box", cex.title = 2.8, rotation.title = 0, showTitle = T, fontcolor.title = "white", fontsize = 6)
        Gviz::displayPars(track.promoter.max.trx) <- list(featureAnnotation = "id")
    } else {
        track.promoter.max.trx <- Gviz::AnnotationTrack(start = str.prom, width = prom.width, chr = chr.gene, strand = strand, hg = hg, name = "Promoter", id = direction.arrow, showTitle = T, col = "mediumblue", fill = "mediumblue", shape = "box", cex.title = 2.8, rotation.title = 0, showTitle = T, fontcolor.title = "white", fontsize = 6)
        Gviz::displayPars(track.promoter.max.trx) <- list(featureAnnotation = "id")
    }


	#track for tracking trx / trxs
	suppressWarnings(
        if(annotation == "ucsc"){
            if(show.all.transcripts){
              if(blackandwhite){
                track.gene.transcript <- Gviz::GeneRegionTrack(txdb, genome = hg, chromosome = chr.gene, name = symbol, start = str.gene, end = end.gene)
                Gviz::displayPars(track.gene.transcript) <- list(transcriptAnnotation = "transcript", fontcolor.group = "black", col = "gray20", just.group = "above", rotation.title = 0, fontcolor.item = "black", col.line = "black", fill = "gray61", cex.title = 0.7, fontsize.group = 25, fontsize = 30)
              }else{
                track.gene.transcript <- Gviz::GeneRegionTrack(txdb, genome = hg, chromosome = chr.gene, name = symbol, start = str.gene, end = end.gene)
                Gviz::displayPars(track.gene.transcript) <- list(transcriptAnnotation = "transcript", fontcolor.group = "black", col = "forestgreen", just.group = "above", rotation.title = 0,  fontcolor.item = "black", col.line = "black", fill = "forestgreen", cex.title = 0.7, fontsize.group = 25, fontsize = 30)
              }
            }else{
              if(blackandwhite){
                track.gene.transcript <- Gviz::GeneRegionTrack(txdb, genome = hg, chromosome = chr.gene, name = symbol, start = str.gene, end = end.gene)
                Gviz::displayPars(track.gene.transcript) <- list(transcriptAnnotation = "trascript", collapseTranscripts = "longest", fontcolor.group = "black", col = "gray20", just.group = "above", rotation.title = 0,  fontcolor.item = "black", col.line = "black", fill = "gray61", cex.title = 0.7, fontsize.group = 25, fontsize = 30)
              }else{
                track.gene.transcript <- Gviz::GeneRegionTrack(txdb, genome = hg, chromosome = chr.gene, name = symbol, start = str.gene, end = end.gene)
                Gviz::displayPars(track.gene.transcript) <- list(transcriptAnnotation = "trascript", collapseTranscripts = "longest", fontcolor.group = "black", col = "forestgreen", rotation.title = 0,  just.group = "above", fontcolor.item = "black", col.line = "black", fill = "forestgreen", cex.title = 0.7, fontsize.group = 25, fontsize = 30)
              }}
            } else {
              b = tab.allTrack.for.Gviz.Ensdb[which(tab.allTrack.for.Gviz.Ensdb$symbol == symbol & tab.allTrack.for.Gviz.Ensdb$seqnames == gsub("chr","",chr.gene)),]
              track.gene.transcript = Gviz::GeneRegionTrack(b, genome = hg, chromosome = chr.gene, name = symbol)
              if(show.all.transcripts){
                if(blackandwhite){
                    Gviz::displayPars(track.gene.transcript) <- list(transcriptAnnotation = "transcript", fontcolor.group = "black", col = "gray20", fill = "gray61", just.group = "above", rotation.title = 0,  cex.title = 0.7, fontsize.group = 22, fontsize = 30)
                }else{
                    Gviz::displayPars(track.gene.transcript) <- list(transcriptAnnotation = "transcript", fontcolor.group = "black", col = "forestgreen", fill = "forestgreen", just.group = "above", rotation.title = 0,  cex.title = 0.7, fontsize.group = 22, fontsize = 30)
                }
            }else{
              if(blackandwhite){
                  track.gene.transcript <- Gviz::GeneRegionTrack(b, genome = hg, chromosome = chr.gene, name = symbol)
                  Gviz::displayPars(track.gene.transcript) <- list(transcriptAnnotation = "transcript", collapseTranscripts = "longest", fontcolor.group = "black", col = "gray20", fill = "gray61", rotation.title = 0,  just.group = "above", cex.title = 0.7, fontsize.group = 22, fontsize = 30)
              }else{ 
                  track.gene.transcript <- Gviz::GeneRegionTrack(b, genome = hg, chromosome = chr.gene, name = symbol)
                  Gviz::displayPars(track.gene.transcript) <- list(transcriptAnnotation = "transcript", collapseTranscripts = "longest", col = "forestgreen", fill = "forestgreen", rotation.title = 0,  fontcolor.group = "black", just.group = "above", cex.title = 0.7, fontsize.group = 22, fontsize = 30)
              }
            }
        }
	)

	#track for plot the lenght of meth segment /segments
	met.reg.width = c()
	group = c()
	for(i in 1: nrow(tab.temp)){
	cur.reg = tab.temp[i,"width"]
	met.reg.width = c(met.reg.width,cur.reg)
	}
	group = paste0(as.character(met.reg.width),sep = " ","bp")
	if(blackandwhite){
        col = "white"
	}else{
        col = "whitesmoke"
	}
	if(is.null(beta2.name)){
        if(length(which(end.start < 5000)) >= 1 & length(which(duplicated(met.reg.width)))>0){
            track.met.reg.width = Gviz::AnnotationTrack(start = str.met, width = tab.temp[,"width"], chromosome = chr.gene, strand = "*", genome = hg, name = name.width, rotation.title = 0, fill = col,id = paste0(as.character(met.reg.width),sep = " ","bp"), showTitle = T)
            Gviz::displayPars(track.met.reg.width) = list(shape = "box", featureAnnotation = "id", fontcolor.feature = "black", col = col, rotation.item = rotation.item, showTitle = T, fontsize = 6, cex.title = 2.5)
        }else if(length(which(duplicated(met.reg.width))) == 0){
            track.met.reg.width = Gviz::AnnotationTrack(start = str.met, width = tab.temp[,"width"], chromosome = chr.gene, strand = "*", genome = hg, name = name.width, rotation.title = 0, fill = col,id = paste0(as.character(met.reg.width),sep = " ","bp"), showTitle = T)
            Gviz::displayPars(track.met.reg.width) = list(shape = "box", featureAnnotation = "id", fontcolor.feature = "black", col = col, rotation.item = rotation.item, showTitle = T, fontsize = 6, cex.title = 2.5)
        } else {
            track.met.reg.width = Gviz::AnnotationTrack(start = str.met, width = tab.temp[,"width"], group = group, chromosome = chr.gene, strand = "*", genome = hg, name = name.width, rotation.title = 0, fill = col, showTitle = T)
            Gviz::displayPars(track.met.reg.width) = list(shape = "box", groupAnnotation = "group", fontcolor.group = "black", col = col, rotation.item = rotation.item, showTitle = T, fontsize = 6, fontsize.group = 20, cex.title = 2.5)
        }
    } else {
        track.met.reg.width = Gviz::AnnotationTrack(start = str.met, width = tab.temp[,"width"], chromosome = chr.gene, strand = strand, genome = hg, name = name.width, rotation.title = 0, fill = col, id = paste0(as.character(met.reg.width),sep = " ","bp"), showTitle = T)
        Gviz::displayPars(track.met.reg.width) = list(shape = "box", featureAnnotation = "id", fontcolor.feature = "black", rotation.item = rotation.item, col = col, showTitle = T, fontsize = 6, cex.title = 2.5)
    }

	#get all bvalue from tab.temp, use these for get the track for meth segment/s (below)
	if(!is.null(beta2.name)){
        if(nrow(tab.temp) == 1){
            bvalue.groups = tab.temp[,c('beta1', 'beta2')]
            bvalue = c()
            for(i in 1:length(bvalue.groups)){
            cur.bvalue =	bvalue.groups[[i]]
            bvalue = c(bvalue,cur.bvalue)}
        } else {
            bvalue.groups = tab.temp[,c('beta1', 'beta2')]
            bvalue = c()
            for(i in 1: nrow(bvalue.groups)){
            cur.bvalue = bvalue.groups[i,]
            bvalue = rbind(bvalue,cur.bvalue)
            }
        }
	} else {
	bvalue = c()
	for(i in 1:nrow(tab.temp)){
		cur.bvalue = tab.temp[,'beta1'][i]
		bvalue = c(bvalue,cur.bvalue)}
	}


	## If there is only one segment to plot
	if(nrow(tab.temp) == 1){
        #if there are two bvalues to plot
        if(!is.null(beta2.name)){
            reg.met.granges = data.frame(seqnames = rep(chr.gene,length(bvalue)), start = c(rep(str.met, length(bvalue)), rep(end.met, length(bvalue))), end=	c(rep(str.met, length(bvalue)), rep(end.met, length(bvalue))))
            for(i in 1:length(bvalue)){
            cur.bvalue.granges = data.frame(bvalue = c(0,rep(bvalue[i],2),0))
            reg.met.granges = cbind(reg.met.granges,cur.bvalue.granges)
            colnames(reg.met.granges)[3+i] = paste0("bvalue",i)
            }
            if(reg.met.granges[,"bvalue1"][2] < reg.met.granges[,"bvalue2"][2]){
                groups =	c(beta2.name,beta1.name)
            } else {
                groups =	c(beta1.name, beta2.name)
            }
            reg.met.granges = GenomicRanges::makeGRangesFromDataFrame(reg.met.granges, keep.extra.columns = T)
            reg.met.track = Gviz::DataTrack(reg.met.granges, groups = groups, genome = hg, name = name1, rotation.title = rotation.title)
            Gviz::displayPars(reg.met.track) = list(ylim = c(0,1.001), yTicksAt = c(0,seq(0.1,0.9,0.1),1), group = groups, col = beta.colors, box.legend = T, cex.axis = cex.axis, cex.legend = 0.7, fontcolor.legend = "black", fontsize = 15 ,cex.title = cex.title.metTrack )
        } else {
            #if there is only one bvalue to plot
            reg.met.granges = data.frame(seqnames = rep(chr.gene,2), start = c(rep(str.met,2), rep(end.met, 2)), end=	c(rep(str.met, 2), rep(end.met, 2)))
            cur.bvalue.granges = data.frame(bvalue = c(0,rep(tab.temp[,'beta1'],c(tab.temp[,'beta1']),2),0))
            reg.met.granges = cbind(reg.met.granges,cur.bvalue.granges)
            reg.met.granges = reg.met.granges[-3,]
            reg.met.granges[2,3] = reg.met.granges[3,2]
            reg.met.granges = GenomicRanges::makeGRangesFromDataFrame(reg.met.granges, keep.extra.columns= TRUE)
            reg.met.track = Gviz::DataTrack(reg.met.granges, genome = hg, name = name2, rotation.title = rotation.title )
            if(length(which(tab.temp$beta1<0)) >0 & length(which(tab.temp$beta1>0)) == 0){
                ylim = c(-1.001,0)
            } else if(length(which(tab.temp$beta1<0)) >0 & length(which(tab.temp$beta1>0)) > 0){
                ylim = c(-1.001,1.001)
            } else {
                ylim = c(0,1.001)
            }
            Gviz::displayPars(reg.met.track) = list(ylim = ylim, yTicksAt = c(0,seq(0.1,0.9,0.1),1), fill.histogram = beta.colors[1], col.histogram = beta.colors[1], box.legend = T ,	cex.axis = cex.axis, cex.legend = 0.7, fontcolor.legend = "black", fontsize = 15, cex.title = cex.title.metTrack)
        }


	    #If path isn't NULL, get the pdf file of the plot (on the path selected)
        if(!is.null(path) & is.null(coord.zoom)){
            pdf(paste0(path,title,".pdf"), height = height.pdf, width = width.pdf)
            p = Gviz::plotTracks(list(track.gene.chromosome,track.gene.genome.axis,track.gene.CpG.island,track.promoter.max.trx,track.gene.transcript,reg.met.track,track.met.reg.width),
                from = from, to = to,
                cex = 2,
                fontface.main = "Times New Roman",
                title.width = title.width.save ,
                type = type,
                background.title = background.title,
                background.panel = background.panel,
                background.legend= NULL
            )
            dev.off()
        } else {
            p <- Gviz::plotTracks(list(track.gene.chromosome, track.gene.genome.axis, track.gene.CpG.island, track.promoter.max.trx, track.gene.transcript, reg.met.track, track.met.reg.width),
                from = from, to = to,
                cex = 2.8,
                fontface.main = "Times New Roman",
                title.width = title.width,
                type = type,
                background.title = background.title,
                background.panel = background.panel,
                background.legend = NULL
            )
        }
	## If there is more than one meth segment to track
	} else {
	#if there are more than one bvalues to plot
	    if(!is.null(beta2.name)){
            reg.met.granges = data.frame()
            for(i in 1 : nrow(tab.temp)){
                cur.reg.met = data.frame(seqnames = rep(chr.gene,4) , start = c(rep(str.met[i],2), rep(end.met[i], 2)), end=	c(rep(str.met[i],2), rep(end.met[i], 2)))
                for(t in 1 : length(bvalue)){
                    cur.bvalue = as.data.frame(c(-1,bvalue[i,][[t]],bvalue[i,][[t]],-1))
                    colnames(cur.bvalue)[1] = paste0("bvalue",t)
                    cur.reg.met = cbind(cur.reg.met,cur.bvalue)
                }
                reg.met.granges = rbind(reg.met.granges,cur.reg.met)
            }
            if(cur.reg.met[,"bvalue1"][2] < cur.reg.met[,"bvalue2"][2]){
                groups =	c(colnames(tab.temp)[grep(beta2.name,colnames(tab.temp))], colnames(tab.temp)[grep(beta1.name,colnames(tab.temp))])
                beta.colors = c(beta.colors[2],beta.colors[1])
            } else {
                groups =	c(colnames(tab.temp)[grep(beta1.name,colnames(tab.temp))], colnames(tab.temp)[grep(beta2.name,colnames(tab.temp))])
            }
            reg.met.granges = GenomicRanges::makeGRangesFromDataFrame(reg.met.granges, keep.extra.columns = T)
            reg.met.track = Gviz::DataTrack(reg.met.granges, groups = groups, genome = hg, name =	name1, rotation.title = rotation.title )
            Gviz::displayPars(reg.met.track) = list(ylim = c(0,1.001),yTicksAt = c(0,seq(0.1,0.9,0.1),1), group = groups, col = beta.colors, box.legend = T, cex.axis = cex.axis, cex.legend =0.7, fontcolor.legend = "black", fontsize = 15, cex.title = cex.title.metTrack)
        } else {
            #if there is only one bvalue to plot
            reg.met.granges = data.frame()
            for(i in 1 : nrow(tab.temp)){
                cur.reg.met = data.frame(seqnames = rep(chr.gene,3) , start = c(rep(str.met[i],3), rep(end.met[i], 3)), end =	c(rep(str.met[i],3), rep(end.met[i], 3)), bvalue = c(0,tab.temp[,'beta1'][i]))
                cur.reg.met = cur.reg.met[-6,]
                cur.reg.met[3,4] = cur.reg.met[2,4]
                cur.reg.met[3,3] = cur.reg.met[4,2]
                cur.reg.met = cur.reg.met[3,]
                reg.met.granges = rbind(reg.met.granges,cur.reg.met)
            }
            reg.met.granges = GenomicRanges::makeGRangesFromDataFrame(reg.met.granges, keep.extra.columns = T)
            reg.met.track = Gviz::DataTrack(reg.met.granges, genome = hg, name = name2, rotation.title = rotation.title )
            if(length(which(tab.temp$beta1<0)) >0 & length(which(tab.temp$beta1>0)) == 0){
                ylim = c(-1.001,0)
            } else if(length(which(tab.temp$beta1<0)) >0 & length(which(tab.temp$beta1>0)) > 0) {
                ylim = c(-1.001,1.001)
            } else {
                ylim = c(0,1.001)
            }
            Gviz::displayPars(reg.met.track) = list(ylim = ylim, yTicksAt = c(0,seq(0.1,0.9,0.1),1), fill.histogram = beta.colors[1], col.histogram = beta.colors[1], box.legend = T, cex.axis = cex.axis, cex.legend =0.7, fontcolor.legend = "black", fontsize = 15,cex.title = cex.title.metTrack)
        }


        # If path isn't NULL, get the pdf file of the plot (on the path selected)
        if(!is.null(path) & is.null(coord.zoom)){
            pdf(paste0(path,title,".pdf"), height = height.pdf, width = width.pdf)
            p = Gviz::plotTracks(list(track.gene.chromosome,track.gene.genome.axis,track.gene.CpG.island,track.promoter.max.trx,track.gene.transcript,reg.met.track,track.met.reg.width),
                        from = from, to = to,
                        cex = 2,
                        fontface.main = "Times New Roman",
                        title.width = title.width.save ,
                        type = type,
                        background.title = background.title,
                        background.panel = background.panel,
                        background.legend= NULL
            )
            dev.off()
        } else {
            p = Gviz::plotTracks(list(track.gene.chromosome, track.gene.genome.axis, track.gene.CpG.island, track.promoter.max.trx, track.gene.transcript, reg.met.track, track.met.reg.width),
                from = from, to = to,
                cex = 2.8,
                fontface.main = "Times New Roman",
                title.width = title.width,
                type = type,
                background.title = background.title,
                background.panel = background.panel,
                background.legend = NULL
            )
        }
	}

	#If coord.zoom isn't NULL, the plot is zoomed to these coordinates
	if(is.null(coord.zoom)){
    	zoom = NULL
	} else {
        #if the zoom is in a plot with only one meth segment
        if(nrow(tab.temp) == 1){
            #if smart zoom is TRUE, adjust the value of coord.zoom ranges
            if(smartzoom){
                from.zoom = coord.zoom[1]
                to.zoom = coord.zoom[2]
                if(to.zoom - from.zoom > 5000){
                    from.zoom = from.zoom - 3000
                    to.zoom = to.zoom + 3000
                }else if(to.zoom - from.zoom < 500){
                    from.zoom = from.zoom -100
                    to.zoom = to.zoom + 100
                } else {
                    from.zoom = from.zoom - 1000
                    to.zoom = to.zoom + 1000
                }
                #if smart zoom is FALSE, the zoom is precisely respect coord.zoom ranges
            } else {
                from.zoom = coord.zoom[1] -1
                to.zoom = coord.zoom[2] +1
            }


            #If path isn't NULL, get the pdf file of the plot (on the path selected)
            if(!is.null(path)){
                pdf(paste0(path,title.zoom, sep = ".", from.zoom, sep = ".", to.zoom,".pdf"), height = height.pdf, width = width.pdf)
                z = Gviz::plotTracks(list(track.gene.chromosome,track.gene.genome.axis,track.gene.CpG.island,track.promoter.max.trx,track.gene.transcript,reg.met.track,track.met.reg.width),
                                    from = from.zoom,
                                    to = to.zoom,
                                    cex = 2,
                                    fontface.main = "Times New Roman",
                                    title.width = title.width.save,
                                    type = type,
                                    background.title = background.title,
                                    background.panel = background.panel,
                                    background.legend= NULL
                )
                dev.off()
            } else {
                z = Gviz::plotTracks(list(track.gene.chromosome, track.gene.genome.axis, track.gene.CpG.island, track.promoter.max.trx, track.gene.transcript, reg.met.track, track.met.reg.width),
                    from = from.zoom, to = to.zoom,
                    cex = 2.8,
                    fontface.main = "Times New Roman",
                    title.width = title.width,
                    type = type,
                    background.title = background.title,
                    background.panel = background.panel,
                    background.legend = NULL
                )
            }

        #if the zoom is in a plot with two meth segments
        } else {
            #if smart zoom is TRUE
            if(smartzoom){
                from.zoom = coord.zoom[1]
                to.zoom = coord.zoom[2]
                if(to.zoom - from.zoom > 5000){
                    from.zoom = from.zoom - 3000
                    to.zoom = to.zoom + 3000
                }else if(to.zoom - from.zoom < 500){
                    from.zoom = from.zoom -100
                    to.zoom = to.zoom + 100
                } else {
                    from.zoom = from.zoom - 1000
                    to.zoom = to.zoom + 1000
                }
                #if smart zoom is FALSE
            } else {
                from.zoom = coord.zoom[1] -1
                to.zoom = coord.zoom[2] +1
            }

            #If path isn't NULL, get the pdf file of the plot (on the path selected)
            if(!is.null(path)){
                pdf(paste0(path,title.zoom,from.zoom, sep = ".", to.zoom,".pdf"), height = height.pdf, width = width.pdf)
                z = Gviz::plotTracks(list(track.gene.chromosome,track.gene.genome.axis,track.gene.CpG.island,track.promoter.max.trx,track.gene.transcript,reg.met.track,track.met.reg.width),
                                from = from.zoom,
                                to = to.zoom,
                                cex = 2,
                                fontface.main = "Times New Roman",
                                title.width = title.width.save,
                                type = type,
                                background.title = background.title,
                                background.panel = background.panel,
                                background.legend= NULL
                )
                dev.off()
            } else {
                z = Gviz::plotTracks(list(track.gene.chromosome, track.gene.genome.axis, track.gene.CpG.island, track.promoter.max.trx, track.gene.transcript, reg.met.track, track.met.reg.width),
                    from = from.zoom, to = to.zoom,
                    cex = 2.8,
                    fontface.main = "Times New Roman",
                    title.width = title.width,
                    type = type,
                    background.title = background.title,
                    background.panel = background.panel,
                    background.legend = NULL
                )
            }
        }
	}
}