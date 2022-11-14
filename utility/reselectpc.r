#!/usr/bin/env Rscript

library(optparse)

ref_reselectpc <- function(pcadir, ref) {

	pcafiles <- unlist(strsplit(list.files(pcadir,pattern=".pc.txt"),".pc.txt"))
	for(j in 1:length(pcafiles)) {
		pcadf <- read.table(paste0(pcadir,"/",pcafiles[j],".pc.txt"), h=T, as.is=T)
		rownames(pcadf) <- paste0(pcadf$chr,"_",pcadf$start)
		ref.chr <- ref[ref$chr==as.character(pcafiles[j]),]
		pcadf[,"ref"] <- ref.chr[rownames(pcadf),"signal"]
		pcavec <- data.frame(val=cor(pcadf[,5:ncol(pcadf)])[1:(ncol(pcadf)-5),(ncol(pcadf)-4)])
		pcavec[,"sign"] <- sign(pcavec$val)
		pcavec$val <- abs(pcavec$val)
		pcavec <- pcavec[order(-pcavec$val),]
		pcavec[,"chr"] <- pcafiles[j]
		pcadf <- pcadf[,c("chr","start","end",rownames(pcavec)[1])]
		pcadf[,4] <- pcadf[,4] * pcavec$sign[1]
		write.table(pcavec, file=paste0(pcadir,"/",pcafiles[j],".pc.ref.details"), col.names=T, row.names=T, sep="\t", quote=F)
		write.table(pcadf, file=paste0(pcadir,"/",pcafiles[j],".pc.bedGraph"), col.names=F, row.names=F, sep="\t", quote=F)
	}
}


man_reselectpc <- function(pcadir, chrom, pcnum, pctype) {

	pcafile  <- list.files(paste0(pcadir,"/",pctype,"_pca"), pattern="_mat", full.names=T)
	cmd <- paste0("cp ",pcafile,"/",chrom,".PC",pcnum,".bedGraph ",pcafile,"/",chrom,".pc.bedGraph")
	cat ("Running ",cmd,"\n")
 	system(cmd, wait=T)
    	cat ("Rewrote ",paste0(pcafile,"/",chrom,".pc.bedGraph"),"\n")
}

#Get the options 
option_list = list(
	
	make_option(c("--reselect"), type="character", help="There are two options to select from\n
		ref : Users should provide a reference bedGraph file with signal and the PC with highest correlation will be selected
		man : Mannual selection of sample wise PC\n
		[ref/man]\n"),

	make_option(c("--rfile"), type="character", default=NA, help="A reference two column text file, e.g.\n
		<_pca folder>\t<reference_bedGraph>\n
		[Must when reselect option is set to ref]\n"),

	make_option(c("--sample"), type="character", default=NA, help="The pca sample folder name [Must when reselect option is set to man]\n"),

	make_option(c("--chr"), type="character", default=NA, help="The chromosome name for which the PC should be reselected [Must when reselect option is set to man]\n"),

	make_option(c("--pc"), type="integer", default=NA, help="The PC number that should be used [Must when reselect option is set to man]\n"),

	make_option(c("--pctype"), type="character", default=NA, help="Reselect intra or inter PC [Must when reselect option is set to man, intra/inter]\n")
)

opt <- parse_args(OptionParser(option_list=option_list))
reselect <- as.character(opt$reselect)
rfile    <- as.character(opt$rfile)
sname    <- as.character(opt$sample) 
chrom    <- as.character(opt$chr)
pcnum    <- as.integer(opt$pc)
pctype   <- as.character(opt$pctype)

if (reselect == "ref" & !is.na(rfile)) {
	reselect.df   <- read.table(rfile, h=F, as.is=T)
	colnames(reselect.df) <- c("dir","ref")
	for(i in 1:nrow(reselect.df)) {
		if (dir.exists(paste0(reselect.df$dir[i],"/intra_pca"))) {
			pcadir <- list.files(paste0(reselect.df$dir[i],"/intra_pca"),pattern="_mat")
			ref.df <- as.data.frame(data.table::fread(reselect.df$ref[i], h=F))
			ref.df <- ref.df[,1:4]
			colnames(ref.df) <- c("chr","start","end","signal")
			rownames(ref.df) <- paste0(ref.df$chr,"_",ref.df$start)
			ref_reselectpc(paste0(reselect.df$dir[i],"/intra_pca/",pcadir), ref.df)
		}
	}
} else if (reselect == "man" & !is.na(sname) & !is.na(chrom) & !is.na(pcnum) & !is.na(pctype)) {
	man_reselectpc(sname, chrom, pcnum, pctype)
}
