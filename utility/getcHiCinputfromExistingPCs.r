#!/usr/bin/env Rscript

library(optparse)

#Get the options
option_list = list(

        make_option(c("--input"), type="character", help="This is a five column text file with the following fields\n
                <pc bedGraph path>\t<pc type>\t<replicate name>\t<sample name>
                
		First column is the path to the bedGraph file with PC information (4th column of the bedGraph file)
		Second column can have either intra or inter
		Thrid column is the replicate name
		Fourth column is the sample name
		
		Note: replicate names and sample names should be different \n")
)
opt <- parse_args(OptionParser(option_list=option_list))

pc_file <- read.table(opt$input, h=F, as.is=T)
colnames(pc_file) <- c("pcfile", "pctype", "repname", "samplename")
head(pc_file)
input_dchic <- list()
cat ("Creating folders for each replicate\n")
for(i in 1:nrow(pc_file)) {
	if (!dir.exists(paste0(pc_file$repname[i],"_HiCPro_pca"))) {
		dir.create(paste0(pc_file$repname[i],"_pca"))
		dir.create(paste0(pc_file$repname[i],"_pca/",pc_file$pctype[i],"_pca"))
		dir.create(paste0(pc_file$repname[i],"_pca/",pc_file$pctype[i],"_pca/",pc_file$repname[i],"_mat"))
	}
	pcbed <- read.table(pc_file$pcfile[i], h=F, as.is=T)
	colnames(pcbed) <- c("chr","start","end","pc")
	print (head(pcbed))
	chroms<- unique(as.character(pcbed$chr))
	print (chroms)
	for(j in 1:length(chroms)) {
		write.table(pcbed[pcbed$chr==chroms[j],], file=paste0(pc_file$repname[i],"_pca/",pc_file$pctype[i],"_pca/",pc_file$repname[i],"_mat/",chroms[j],".pc.bedGraph"), row.names=F, col.names=F, quote=F, sep="\t")
	} 
  	input_dchic[[i]] <- data.frame(mat=paste0(pc_file$repname,".mat"),
				 bed=paste0(pc_file$repname,".bed"),
				 rep=pc_file$repname,
			    	 sam=pc_file$samplename)
}
input_dchic <- unique(do.call(rbind, input_dchic))
write.table(input_dchic, file="input_dcHiC.txt", row.names=F, col.names=F, sep="\t", quote=F)
