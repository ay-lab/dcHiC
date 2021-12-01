#libraries required
library(bigstatsr)
library(Rcpp)
#library(data.table)
#library(parallel)
#library(hashmap)
#library(functionsdchic)
#library(bench)
#library(limma)
#library(depmixS4)
library(optparse)

########################### PCA calculation ###########################
#Generate FBM matrix
mat2fbm <- function(x_list, bin) {
	mat.fbm <- bigstatsr::FBM(bin, bin)
	for(j in 1:(length(x_list)-1)) {
		r.start <- x_list[[j]]$row.start + 1
		c.start <- x_list[[j]]$col.start + 1
		r.end 	<- x_list[[j]]$row.end + 1
		c.end 	<- x_list[[j]]$col.end + 1
		x_list[[j]]$mat[is.na(x_list[[j]]$mat)] <- 0
		if (r.start == c.start) {
			diag(x_list[[j]]$mat) <- 1
			mat.fbm[r.start:r.end, c.start:c.end] <- round(x_list[[j]]$mat,5)
		} else {
			mat.fbm[r.start:r.end, c.start:c.end] <- round(x_list[[j]]$mat,5)
			mat.fbm[c.start:c.end, r.start:r.end] <- round(t(x_list[[j]]$mat),5)
		}
	}
	return(mat.fbm)
}

#Write inter chromosomal interaction to file
extractTrans <- function(l, chr, matpath, bedpath, outpath) {
	functionsdchic::createtransijk(matpath, bedpath, paste0(outpath,"/",chr[l],".inter.txt"), as.character(chr[l]))	
}

#Generate trans matrix and run the pca
ijk2matfunc_trans <- function(l, chr, beddf, path, pc_k = 3, pc_ncores, rmerge = 10000) {

	rows<- nrow(beddf[beddf$chr != as.character(chr[l]),])
	cols<- nrow(beddf[beddf$chr == as.character(chr[l]),])

	#Divide the rows
	bin <- rows
   	bindiv <- round(seq(1, bin, length.out=rmerge))
   	rbins <- c()
   	rbine <- c()
   	rbins[1] <- bindiv[1]
   	rbine[1] <- bindiv[2]
   	k <- 2
	while (k < length(bindiv)) {
		rbins[k] <- bindiv[k] + 1
		rbine[k] <- bindiv[(k+1)] 
		k <- k + 1
	}

	#Divide the colunms
	bin <- cols
   	bindiv <- round(seq(1, bin, length.out=10))
   	cbins <- c()
   	cbine <- c()
   	cbins[1] <- bindiv[1]
   	cbine[1] <- bindiv[2]
   	k <- 2
	while (k < length(bindiv)) {
		cbins[k] <- bindiv[k] + 1
		cbine[k] <- bindiv[(k+1)] 
		k <- k + 1
	}


    #Remove single instance
    trans_mat <- as.data.frame(data.table::fread(paste0(path,"/",chr[l],".inter.txt"), h=F))
    trans_cov <- aggregate(V3 ~ V2, sum, data=trans_mat)
    trans_cov <- trans_cov[trans_cov$V3 > 1,]
    trans_mat <- trans_mat[trans_mat$V2 %in% c(trans_cov$V2),]
    write.table(trans_mat, file=paste0(path,"/",chr[l],".inter.txt"), row.names=F, col.names=F, sep="\t", quote=F)

	#Run the Z transformation and then correlation calculation
	x <- functionsdchic::transmat(paste0(path,"/",chr[l],".inter.txt"), (rbins-1), (rbine-1), (cbins-1), (cbine-1), rows, cols)

	#info <- bench::mark(min_time = Inf, max_iterations=1, time_unit='s', x <- functionsdchic::transmat(paste0(path,"/",chr[l],".inter.txt"), (rbins-1), (rbine-1), (cbins-1), (cbine-1), rows, cols))
	#sink(paste0(path,"/",chr[l],".correlation.system.txt"))
   	#print (cbind(time.seconds=info$min,mem.gb=info$mem_alloc/8e9))
   	#sink()	

	#Covert the matrix into a FBM object
	mat.fbm <- mat2fbm(x$zcor, bin)
	remove  <- which(x$coverage < 10)
	rm(x)

	#Reformat the bedfile if blacklisted regions are provided

	if (ncol(beddf) == 4) {
		bldf <- beddf[beddf$chr == as.character(chr[l]),]
		colnames(bldf) <- c("chr", "start", "end", "index")
		bldf[,"index"] <- 1:nrow(bldf)
		bldf[,"blacklist"] <- 0
   		bldf[remove,"blacklist"] <- 1
		blvec <- as.integer(bldf[bldf$blacklist == 1,]$index)
		mat.fbm <- bigstatsr::as_FBM(mat.fbm[-c(blvec),-c(blvec)])

	} else if (ncol(beddf) == 5) {
		bldf <- beddf[beddf$chr == as.character(chr[l]),]
		colnames(bldf) <- c("chr", "start", "end", "index","blacklist")
		bldf[,"index"] <- 1:nrow(bldf)
   		bldf[remove,"blacklist"] <- 1
		blvec <- as.integer(bldf[bldf$blacklist == 1,]$index)
		mat.fbm <- bigstatsr::as_FBM(mat.fbm[-c(blvec),-c(blvec)])
	}


	if (file.exists(paste0(path,"/",chr[l],".precmat.txt"))) {
		system(paste0("rm ",path,"/",chr[l],".precmat.txt"), wait=T)
	}
	if (file.exists(paste0(path,"/",chr[l],".cmat.txt"))) {
		system(paste0("rm ",path,"/",chr[l],".cmat.txt"), wait=T)
	}
	bigstatsr::big_write(mat.fbm,  every_nrow = 100, progress = interactive(), file=paste0(path,"/",chr[l],".precmat.txt"))

	bin <- nrow(mat.fbm)
   	bindiv <- round(seq(1, bin, length.out=10))
   	start  <- c()
   	end	   <- c()
   	start[1] <- bindiv[1]
   	end[1]   <- bindiv[2]
   	k <- 2
	while (k < length(bindiv)) {
		start[k] <- bindiv[k] + 1
		end[k]   <- bindiv[(k+1)] 
		k <- k + 1
	}

	#Calculate correlation matrix again on the previously calculated correlation matrix
	#The last 0 value will not return the total coverage of the bins, 1 will. 
	x <- functionsdchic::oe2cor(functionsdchic::fbm2mat(mat.fbm), (start-1), (end-1), 1, 0)

	#Convert it to FMB object
	mat.fbm <- mat2fbm(x, nrow(mat.fbm))	

	#Perform SVD on the correlation matrix
   	obj.svd <- bigstatsr::big_randomSVD(mat.fbm, fun.scaling = bigstatsr::big_scale(center = FALSE, scale = FALSE), k = pc_k, ncores = pc_ncores)
	#info <- bench::mark(min_time = Inf, max_iterations=1, time_unit='s', obj.svd <- bigstatsr::big_randomSVD(mat.fbm, fun.scaling = bigstatsr::big_scale(center = FALSE, scale = FALSE), k = pc_k, ncores = pc_ncores))
	#sink(paste0(path,"/",chr[l],".svd.system.txt"))
   	#print (cbind(time.seconds=info$min,mem.gb=info$mem_alloc/8e9))
   	#sink()

   	# Get the PCs (R-mode)
	pc <- x[[length(x)]]$zmat %*% obj.svd$v
	bldf <- bldf[bldf$blacklist == 0,]
	bldf <- cbind(bldf[,1:4],pc)
	colnames(bldf) <- c("chr","start","end","index",paste0("PC",1:pc_k))
   	
   	write.table(bldf, file=paste0(path,"/",chr[l],".pc.txt"), row.names=F, sep="\t", quote=F)
	saveRDS(obj.svd, file=paste0(path,"/",chr[l],".svd.rds"))
	bigstatsr::big_write(mat.fbm,  every_nrow = 100, progress = interactive(), file=paste0(path,"/",chr[l],".cmat.txt"))
	rm(x)
	rm(obj.svd)
}

#Read the files and create trans interation files
readfilesinter <- function(i, df, pcout, cthrd, pthrd, rmergeno = 10000, diroverwrite=FALSE) {
	
	if (!dir.exists(paste0(df$prefix[i],"_pca"))) {
		dir.create(paste0(df$prefix[i],"_pca"))
	}
	if (!dir.exists(paste0(df$prefix[i],"_pca/","inter_pca"))) {
		dir.create(paste0(df$prefix[i],"_pca/","inter_pca"))
	}
	sink(file=paste0(df$prefix[i],"_pca/",df$prefix[i],"_pca.inter.log"))
	cat ("Generating chromosome-wise trans interactions files\n")
	if (!dir.exists(paste0(df$prefix[i],"_pca/","inter_pca/",df$prefix[i],"_mat")) | diroverwrite == TRUE) {
		dir.create(paste0(df$prefix[i],"_pca/","inter_pca/",df$prefix[i],"_mat"))
		bed   <- read.table(normalizePath(df$bed[i]), h=F, as.is=T)
		if (ncol(bed) == 4) {
			colnames(bed) <- c("chr","start","end","index")
		} else if (ncol(bed) == 5) {
			colnames(bed) <- c("chr","start","end","index","blacklist")
		}
		chrom <- unique(as.character(bed[,1]))
		cat ("Found ",chrom," unique chromosomes\n")
		if (cthread > 1) {
			cl_cthread <- parallel::makeCluster(cthrd)
			parallel::clusterExport(cl_cthread, c("extractTrans","ijk2matfunc_trans","mat2fbm"))
			parallel::parLapply(cl_cthread, 1:length(chrom), extractTrans, chrom, normalizePath(df$mat[i]), normalizePath(df$bed[i]), paste0(df$prefix[i],"_pca/","inter_pca/",df$prefix[i],"_mat"))
			parallel::parLapply(cl_cthread, 1:length(chrom), ijk2matfunc_trans, chrom, bed, paste0(df$prefix[i],"_pca/","inter_pca/",df$prefix[i],"_mat"), pcout, pthrd, rmergeno)
			parallel::stopCluster(cl_cthread)
		} else {
			lapply(1:length(chrom), extractTrans, chrom, normalizePath(df$mat[i]), normalizePath(df$bed[i]), paste0(df$prefix[i],"_pca/","inter_pca/",df$prefix[i],"_mat"))
			lapply(1:length(chrom), ijk2matfunc_trans, chrom, bed, paste0(df$prefix[i],"_pca/","inter_pca/",df$prefix[i],"_mat"), pcout, pthrd, rmergeno)
		}
	} else {
		cat (paste0(df$prefix[i],"_pca/","inter_pca/",df$prefix[i],"_mat")," folder exists!\n")
	}
	sink()
}


# Total interaction
expectedInteraction <- function(i, dist_param, bed_nrow, start_resolution) {
 
  	vecInt <- bed_nrow$nrow - (dist_param$dist[i]/start_resolution)
  	vecInt <- sum(vecInt[vecInt > 0])
  	return(vecInt)
}

## This function converts the Hi-C txt file to matrix file
ijk2matfunc_cis <- function(i, chr, path, pc_k = 3, pc_ncores) {

   	df  <- data.table::fread(paste0(path,"/",chr[i],".txt"), h=T)
   	bed <- read.table(paste0(path,"/",chr[i],".bed"), h=T, as.is=T)
   	mat <- functionsdchic::ijk2mat(as.matrix(df[,c("A","B","WeightOE")]), nrow(bed), nrow(bed))
   	if (ncol(bed) == 4) {
		colnames(bed) <- c("chr","start","end","index")
		remove <- which(apply(mat, 1, sum) < 3)
   		bed[,"blacklist"] <- 0
   		bed[remove,"blacklist"] <- 1
	} else if (ncol(bed) == 5) {
		colnames(bed) <- c("chr","start","end","index","blacklist")
		remove <- which(apply(mat, 1, sum) < 3)
   		bed[,"blacklist"] <- 0
   		bed[remove,"blacklist"] <- 1
	}
   	bin <- nrow(bed)
   	bindiv <- round(seq(1, bin, length.out=10))
   	start  <- c()
   	end	   <- c()
   	start[1] <- bindiv[1]
   	end[1]   <- bindiv[2]
   	k <- 2
	while (k < length(bindiv)) {
		start[k] <- bindiv[k] + 1
		end[k]   <- bindiv[(k+1)] 
		k <- k + 1
	}
   	
    
    #The 1 at the instructs to send the zscore matrix, while 0 will not
    #The last 0 value will not return the total coverage of the bins, 1 will. 
   	x <- functionsdchic::oe2cor(mat, (start-1), (end-1), 1, 0)
    rm(mat)

	#sink(paste0(path,"/",chr[i],".correlation.system.txt"))
   	#print (cbind(time.seconds=info$min,mem.gb=info$mem_alloc/8e9))
   	#sink()	

   	#Covert the matrix into a FBM object
	mat.fbm <- mat2fbm(x, bin)
	rm(x)

	if (ncol(bed) == 5) {
		bldf <- bed[bed$chr == as.character(chr[i]),]
		bldf[,"index"] <- 1:nrow(bldf)
		blvec <- as.integer(bldf[bldf$blacklist == 1,]$index)
		if (nrow(bldf[bldf$blacklist == 1,]) > 0) {
			mat.fbm <- bigstatsr::as_FBM(mat.fbm[-c(blvec),-c(blvec)])
		}
	}

	if (file.exists(paste0(path,"/",chr[i],".precmat.txt"))) {
		system(paste0("rm ",path,"/",chr[i],".precmat.txt"), wait=T)
	}
	if (file.exists(paste0(path,"/",chr[i],".cmat.txt"))) {
		system(paste0("rm ",path,"/",chr[i],".cmat.txt"), wait=T)
	}
	bigstatsr::big_write(mat.fbm,  every_nrow = 100, progress = interactive(), file=paste0(path,"/",chr[i],".precmat.txt"))

	bin <- nrow(mat.fbm)
   	bindiv <- round(seq(1, bin, length.out=10))
   	start  <- c()
   	end	   <- c()
   	start[1] <- bindiv[1]
   	end[1]   <- bindiv[2]
   	k <- 2
	while (k < length(bindiv)) {
		start[k] <- bindiv[k] + 1
		end[k]   <- bindiv[(k+1)] 
		k <- k + 1
	}

	#Calculate correlation matrix again on the previously calculated correlation matrix
	#The last 0 value will not return the total coverage of the bins, 1 will. 
	x <- functionsdchic::oe2cor(functionsdchic::fbm2mat(mat.fbm), (start-1), (end-1), 1, 0)

	#Convert it to FMB object
	mat.fbm <- mat2fbm(x, nrow(mat.fbm))

   	#info <- bench::mark(min_time = Inf, max_iterations=1, time_unit='s', obj.svd <- bigstatsr::big_randomSVD(mat.fbm, fun.scaling = bigstatsr::big_scale(center = FALSE, scale = FALSE), k = pc_k, ncores = pc_ncores))
	#sink(paste0(path,"/",chr[i],".svd.system.txt"))
   	#print (cbind(time.seconds=info$min,mem.gb=info$mem_alloc/8e9))
   	#sink()

   	#Perform SVD on the correlation matrix
   	obj.svd <- bigstatsr::big_randomSVD(mat.fbm, fun.scaling = bigstatsr::big_scale(center = FALSE, scale = FALSE), k = pc_k, ncores = pc_ncores)

	# Get the PCs
	pc <- x[[length(x)]]$zmat %*% obj.svd$v
	
	#Reformat the bed file 
   	if (ncol(bed) == 4) {
   		bedchr <- cbind(bed[bed$chr == as.character(chr[i]),1:4],pc)
   		colnames(bedchr) <- c("chr","start","end","index",paste0("PC",1:pc_k))	
   	} else if (ncol(bed) == 5) {
   		bedchr = cbind(bldf[bldf$blacklist == 0,1:4],pc)
   		print (head(bedchr))
   		colnames(bedchr) <- c("chr","start","end","index",paste0("PC",1:pc_k))
   	}
   	write.table(bedchr, file=paste0(path,"/",chr[i],".pc.txt"), row.names=F, sep="\t", quote=F)
	saveRDS(obj.svd, file=paste0(path,"/",chr[i],".svd.rds"))
	bigstatsr::big_write(mat.fbm,  every_nrow = 100, progress = interactive(), file=paste0(path,"/",chr[i],".cmat.txt"))
	rm(x)
	rm(obj.svd)
}

#Read the files and create cis O/E files
readfilesintra <- function(i, df, pcout, ebackground, cthrd, pthrd, diroverwrite=FALSE) {

	mat <- normalizePath(df$mat[i])
	bed <- normalizePath(df$bed[i])

	if (!dir.exists(paste0(df$prefix[i],"_pca"))) {
		dir.create(paste0(df$prefix[i],"_pca"))
	}
	if (!dir.exists(paste0(df$prefix[i],"_pca/","intra_pca"))) {
		dir.create(paste0(df$prefix[i],"_pca/","intra_pca"))
	}

	#sink(file=paste0(df$prefix[i],"_pca/",df$prefix[i],"_pca.intra.log"))
	if (!dir.exists(paste0(df$prefix[i],"_pca/","intra_pca/",df$prefix[i],"_mat")) | diroverwrite == TRUE) {
		cat ("Reading ",mat," Hi-C matrix file\n")
  		cat ("Reading ",bed," Hi-C bed file\n")
 		mat <- data.table::fread(mat, h=F)
 		bed <- read.table(bed, h=F, as.is=T)
 		colnames(mat) <- c("A","B","Weight")
 		if (ncol(bed) == 4) {
			colnames(bed) <- c("chr","start","end","index")
		} else if (ncol(bed) == 5) {
			colnames(bed) <- c("chr","start","end","index","blacklist")
		}
		chr.hash <- hashmap::hashmap(bed$index, as.vector(bed$chr))
 		pos.hash <- hashmap::hashmap(bed$index, c(bed$start))
 		mat[,"chr1"] <- chr.hash[[mat$A]]
 		mat[,"pos1"] <- pos.hash[[mat$A]]
 		mat[,"chr2"] <- chr.hash[[mat$B]]
 		mat[,"pos2"] <- pos.hash[[mat$B]]
 		mat <- mat[mat$chr1==mat$chr2,]
 		mat[,"dist"] <- abs(mat$pos1 - mat$pos2) 
 		print (head(mat))

 		if (!file.exists(paste0(df$prefix[i],"_pca/","intra_pca/",df$prefix[i],"_mat"))) {
 			dir.create(paste0(df$prefix[i],"_pca/","intra_pca/",df$prefix[i],"_mat"))
 		}

		if (ebackground == 1) {

			cat ("Calculating expected counts from overall genomic background\n")
 			dist_param <- aggregate(Weight ~ dist, sum, data=mat)
 			print (head(dist_param))  
  		
 			bed_nrow <- bed
 			bed_nrow[,"nrow"] <- 1
 			bed_nrow <- aggregate(nrow ~ chr, sum, data=bed_nrow)
 			expInt <- lapply(c(1:nrow(dist_param)), expectedInteraction, dist_param, bed_nrow, (bed$end[1] - bed$start[1]))  
 			expInt <- do.call(rbind, expInt)
 			dist_param[,"totalpairinteraction"] <- expInt
 			dist_param[,"expcc"] <- dist_param$Weight/dist_param$totalpairinteraction
 			data.table::fwrite(dist_param, file=paste0(df$prefix[i],"_pca/","intra_pca/",df$prefix[i],"_mat/Genome.distparam"), row.names=F, col.names=T, sep="\t", quote=F)

 			expcc_hash <- hashmap::hashmap(dist_param$dist, dist_param$expcc)
 			mat[,"WeightOE"] <- mat$Weight/expcc_hash[[mat$dist]]
 			print (head(mat))
 			print (nrow(mat))
 			print (dist_param[dist_param$expcc < 3,]$dist[1])
 			mat <- mat[mat$dist < dist_param[dist_param$expcc < 3,]$dist[1],]
 			print (nrow(mat))

  			chr.uniq <- unique(bed$chr)
  			for(j in 1:length(chr.uniq)) {
  	 		 	cat ("Writing ",as.character(chr.uniq[j]),".txt file\n")
  	 			b <- bed[bed$chr==chr.uniq[j],]
  	 			b[,"index"] <- 1:nrow(b)
  				index.hash  <- hashmap::hashmap(b$start, b$index)
  				m <- mat[(mat$chr1==chr.uniq[j] & mat$chr2==chr.uniq[j]), c(4,5,6,7,9,1,2,3)]
  				m[,"A"] <- index.hash[[m$pos1]]
  				m[,"B"] <- index.hash[[m$pos2]] 
  		  		data.table::fwrite(m, file=paste0(df$prefix[i],"_pca/","intra_pca/",df$prefix[i],"_mat/",chr.uniq[j],".txt"), row.names=F, col.names=T, sep="\t", quote=F)
  		  	 	data.table::fwrite(b, file=paste0(df$prefix[i],"_pca/","intra_pca/",df$prefix[i],"_mat/",chr.uniq[j],".bed"), row.names=F, col.names=T, sep="\t", quote=F)
  			}
  		} else if (ebackground == 2) {

  			chr.uniq <- unique(bed$chr)
  			chrom_exist <- unlist(strsplit(list.files(paste0(df$prefix[i],"_pca/","intra_pca/",df$prefix[i],"_mat"), pattern="pc.txt"),".pc.txt"))
  			chr.uniq <- setdiff(chr.uniq, chrom_exist)
  			if (length(chr.uniq) > 0) {
  				for(j in 1:length(chr.uniq)) {

  					cat ("Calculating expected counts from chromosome wise background\n")
  					mat_chr <- mat[mat$chr1 == as.character(chr.uniq[j]) & mat$chr2 == as.character(chr.uniq[j]),]
  					bed_nrow <- bed[bed$chr == as.character(chr.uniq[j]),]
  					dist_param <- aggregate(Weight ~ dist, sum, data=mat_chr)
 					print (head(dist_param))
 				
 					bed_nrow[,"nrow"] <- 1
 					bed_nrow <- aggregate(nrow ~ chr, sum, data=bed_nrow)
 					expInt <- lapply(c(1:nrow(dist_param)), expectedInteraction, dist_param, bed_nrow, (bed$end[1] - bed$start[1]))  
 					expInt <- do.call(rbind, expInt)
 					dist_param[,"totalpairinteraction"] <- expInt
 					dist_param[,"expcc"] <- dist_param$Weight/dist_param$totalpairinteraction
 					data.table::fwrite(dist_param, file=paste0(df$prefix[i],"_pca/","intra_pca/",df$prefix[i],"_mat/",chr.uniq[j],".distparam"), row.names=F, col.names=T, sep="\t", quote=F)

 					expcc_hash <- hashmap::hashmap(dist_param$dist, dist_param$expcc)
 					mat_chr[,"WeightOE"] <- mat_chr$Weight/expcc_hash[[mat_chr$dist]]
 					print (head(mat_chr))
 					print (nrow(mat_chr))
 					print (dist_param[dist_param$expcc < 1,]$dist[1])
 					#mat_chr <- mat_chr[mat_chr$dist <= dist_param[dist_param$expcc <= 1,]$dist[1],]
 					#mat_chr <- mat_chr[mat_chr$WeightOE >= 0.1 & mat_chr$WeightOE <= 2,]
 					print (nrow(mat_chr))

 					cat ("Writing ",as.character(chr.uniq[j]),".txt file\n")
 					b <- bed[bed$chr==chr.uniq[j],]
  	 				b[,"index"] <- 1:nrow(b)
  					index.hash  <- hashmap::hashmap(b$start, b$index)
  					m <- mat_chr[, c(4,5,6,7,9,1,2,3)]
  					m[,"A"] <- index.hash[[m$pos1]]
  					m[,"B"] <- index.hash[[m$pos2]] 
  		  			data.table::fwrite(m, file=paste0(df$prefix[i],"_pca/","intra_pca/",df$prefix[i],"_mat/",chr.uniq[j],".txt"), row.names=F, col.names=T, sep="\t", quote=F)
  		  	 		data.table::fwrite(b, file=paste0(df$prefix[i],"_pca/","intra_pca/",df$prefix[i],"_mat/",chr.uniq[j],".bed"), row.names=F, col.names=T, sep="\t", quote=F)
  		  	 		rm(mat_chr)
  		  	 		rm(m)
  		  	 	}
  			}
  		}
  		rm(mat)

  		chrom <- list.files(paste0(df$prefix[i],"_pca/","intra_pca/",df$prefix[i],"_mat"), pattern="*.bed")
  		chrom <- chrom[!grepl("pc",chrom)]
  		chrom <- chrom[!grepl("bedGraph",chrom)]
  		chrom <- gsub(".bed", "", chrom)
  		chrom_exist <- unlist(strsplit(list.files(paste0(df$prefix[i],"_pca/","intra_pca/",df$prefix[i],"_mat"), pattern="pc.txt"),".pc.txt"))
  		print (chrom_exist)
  		print (chrom)
  		cat ("Found ",chrom_exist," files\n")
  		chrom <- setdiff(chrom, chrom_exist)
 		if (length(chrom) > 0) {
  			if (cthread > 1) {
  				cl_cthread <- parallel::makeCluster(cthrd)
				parallel::clusterExport(cl_cthread, c("ijk2matfunc_cis","mat2fbm"))
				parallel::parLapply(cl_cthread, 1:length(chrom), ijk2matfunc_cis, chrom, paste0(df$prefix[i],"_pca/","intra_pca/",df$prefix[i],"_mat"), pcout, pthrd)
				parallel::stopCluster(cl_cthread)
			} else {
				lapply(1:length(chrom), ijk2matfunc_cis, chrom, paste0(df$prefix[i],"_pca/","intra_pca/",df$prefix[i],"_mat"), pcout, pthrd)
			}
			} else {
				cat ("There exists PCA files (pc.txt) inside ",paste0(df$prefix[i],"_pca/","intra_pca/",df$prefix[i],"_mat")," folder. Delete those in order to rerun!\n")
			}
  	}
	#sink()
}
### END ###

########################### PC selection and flipping ###########################

## Find bedtools
.findExecutable <- function(exe, interactive=TRUE) {
	path <- Sys.which(exe)
  	if(all(path=="")) {
    	if(interactive) stop("Executable for ", paste(exe, collapse=" or "), " not found! Please make sure that the software is correctly installed and, if necessary, path variables are set.", call.=FALSE)
    	return(character(0))
  	}
  	path[which(path!="")[1]]
}

## Download the data
datadownload  <- function(genome) {

	if (!file.exists(paste0("cytoBand.txt.gz"))) {
    	cytoband <- paste0("curl -O http://hgdownload.cse.ucsc.edu/goldenPath/",genome,"/database/cytoBand.txt.gz")
    	cat ("Running ",cytoband,"\n")
    	system(cytoband, wait=T)
  	}

  	if (!file.exists(paste0(genome,".chrom.sizes"))) {
  	  genomesize <- paste0("curl -O http://hgdownload.cse.ucsc.edu/goldenPath/",genome,"/bigZips/",genome,".chrom.sizes")
  	  cat ("Running ",genomesize,"\n")
  	  system(genomesize, wait=T)
  	}

  	if (!file.exists(paste0(genome,".fa"))) {
  	  chromfa <- paste0("curl -O http://hgdownload.cse.ucsc.edu/goldenPath/",genome,"/bigZips/",genome,".fa.gz\n")
  	  cat ("Running ",chromfa,"\n")
  	  system(chromfa, wait=T)
  	}

  	if (!file.exists(paste0(genome,".refGene.gtf.gz"))) {
  	  gene <- paste0("curl -O http://hgdownload.cse.ucsc.edu/goldenPath/",genome,"/bigZips/genes/",genome,".refGene.gtf.gz")
  	  cat ("Running ",gene,"\n")
  	  system(gene, wait=T)
  	}
}

## PC selection and flipping core function
pcselectioncore <- function(chr, sam, goldenpath, pc_k, sname, pc_type) {
	d <- 1
	data <- list(clus=list(),cor=list(),chrom=list())
	vals <- list()
	for(i in 1:length(chr)) {
		gcc.values <- list()
		tss.values <- list()
		len.values <- list()
		chrom.list <- list()
		count.vect <- list()
		for(j in 1:length(sam)) {
			cat ("Running ",pc_type," ",chr[i]," in ",sam[j]," sample\n")
			pca <- read.table(paste0(sam[j],"_pca/",pc_type,"_pca/",sam[j],"_mat/",chr[i],".pc.txt"), h=T)
			pca <- pca[,1:(4+pc_k)]
			gcc <- read.table(goldenpath, h=F)
			colnames(gcc) <- c("chr","start","end","gcc","tss")
			rownames(gcc) <- paste0(gcc$chr,"_",gcc$start)
			rownames(pca) <- paste0(pca$chr,"_",pca$start)
			pca[,c("gcc","tss")] <- gcc[rownames(pca),c("gcc","tss")]
			for(k in 1:(ncol(pca)-6)) {
				pca[,c(k+4)] <- sign(cor(pca[,c(k+4)], pca$gcc)) * pca[,c(k+4)]
				write.table(pca[,c(1:3,(k+4))], file=paste0(sam[j],"_pca/",pc_type,"_pca/",sam[j],"_mat/",chr[i],".PC",k,".bedGraph"), row.names=F, col.names=F, sep="\t", quote=F)
			}
			chrom.list[[j]] <- data.frame(pca[,c(5:(ncol(pca)-2))])
			colnames(chrom.list[[j]]) <- paste0(sam[j],".PC",c(1:ncol(chrom.list[[j]])))
			count.vect[[j]] <- rownames(chrom.list[[j]])
			gcc.values[[j]] <- cbind(gcc.cor=round(cor(cbind(chrom.list[[j]],gcc=pca$gcc))[1:ncol(chrom.list[[j]]),(ncol(chrom.list[[j]])+1)],4))
			tss.values[[j]] <- cbind(tss.cor=round(cor(cbind(chrom.list[[j]],tss=pca$tss))[1:ncol(chrom.list[[j]]),(ncol(chrom.list[[j]])+1)],4))
			len.values[[j]] <- cbind(len.cor=round(cor(cbind(chrom.list[[j]],len=1:nrow(pca)))[1:ncol(chrom.list[[j]]),(ncol(chrom.list[[j]])+1)],4))
		}
		gcc.values <- as.data.frame(do.call(rbind, gcc.values))
		tss.values <- as.data.frame(do.call(rbind, tss.values))
		len.values <- as.data.frame(do.call(rbind, len.values))
		gcc.values[,"name"] <- rownames(gcc.values)
		tss.values[,"name"] <- rownames(tss.values)
		len.values[,"name"] <- rownames(len.values)
		vals[[i]]  <- cbind(
				gcc.cor=gcc.values[,1],
				tss.cor=tss.values[,1],
				len.cor=len.values[,1],
				chr=chr[i],
				sample=do.call(rbind,strsplit(gcc.values$name,"[.]"))[,1], 
				pc=do.call(rbind,strsplit(gcc.values$name,"[.]"))[,2])
	
		count.vect <- unlist(count.vect)
		count.vect <- data.frame(number=table(count.vect))
		count.vect <- count.vect[count.vect$number.Freq == length(sam),]

		pc.mat <- list()
		for(j in 1:length(sam)) {
			pc.mat[[j]] <- chrom.list[[j]][as.character(count.vect$number.count.vect),c(1:ncol(chrom.list[[j]]))]
		}
		pc.mat <- as.matrix(do.call(cbind, pc.mat))
		hcl <- hclust(as.dist(round(1-cor(pc.mat),4)))
		cl<- list()
		k <- 1
		h <- 0.05
		while (h < 1) {
			g <- data.frame(group=cutree(hcl, h=h))	
			g[,"name"]  <- rownames(g)
			g[,"name2"] <- paste0(g$group,"-",do.call(rbind,strsplit(g$name,"[.]"))[,1])
			v <- data.frame(member=table(g$name2))
			v$member.Var1 <- do.call(rbind,strsplit(as.character(v$member.Var1),"-"))[,1]
			v <- data.frame(member=table(v$member.Var1))
			s <- v[v$member.Freq == length(sam),]$member.Var1
			if (length(s) > 0) {
				for(u in 1:length(s)) {
					if (length(unique(g[g$group==s[u],]$name)) == length(sam)) {
						cl[[k]] <- g[g$group==s[u],]
						k <- k + 1
					}
				}
			}
			h <- h + 0.025
		}

		if (length(cl) > 0) {
			cl <- as.data.frame(do.call(rbind, cl))
			cl <- unique(cl)
			cl[,"chr"] <- chr[i]
			cl[,"gcc.cor"] <- 0
			cl[,"tss.cor"] <- 0
			for(v in 1:nrow(gcc.values)) {
				cl[cl$name == gcc.values$name[v],"gcc.cor"] <- gcc.values$gcc.cor[v]
				cl[cl$name == gcc.values$name[v],"tss.cor"] <- tss.values$tss.cor[v]
				cl[cl$name == gcc.values$name[v],"len.cor"] <- abs(len.values$len.cor[v])
			}
			data$clus[[d]] <- cl
			data$cor[[d]] <- as.data.frame(aggregate(gcc.cor ~ group + chr, mean, data=cl))
			data$cor[[d]][,"tss.cor"] <- aggregate(tss.cor ~ group + chr, mean, data=cl)$tss.cor
			data$cor[[d]][,"len.cor"] <- aggregate(len.cor ~ group + chr, mean, data=cl)$len.cor
			data$cor[[d]][,"score"] <- apply(data$cor[[d]][,c("tss.cor","gcc.cor")],1,sum)
			data$cor[[d]]$gcc.cor <- round(data$cor[[d]]$gcc.cor, 4)
			data$cor[[d]]$tss.cor <- round(data$cor[[d]]$tss.cor, 4)
			data$cor[[d]]$len.cor <- round(data$cor[[d]]$len.cor, 4)
			data$cor[[d]]$score   <- round(data$cor[[d]]$score, 4)
			data$chrom[[d]] <- chr[i]
			d <- d + 1
		}
	}

	clus.df <- as.data.frame(unique(do.call(rbind,data$clus)))
	cor.df  <- as.data.frame(do.call(rbind,data$cor))
	chr.vec <- unlist(data$chrom)
	vals    <- as.data.frame(do.call(rbind,vals), stringsAsFactors = FALSE)
	vals$gcc.cor <- as.numeric(vals$gcc.cor)
	vals$tss.cor <- as.numeric(vals$tss.cor)
	vals$len.cor <- abs(as.numeric(vals$len.cor))

	write.table(clus.df, file=paste0(sname,"_clus.txt"), row.names=F, col.names=T, sep="\t", quote=F)
	write.table(cor.df, file=paste0(sname,"_cor.txt"), row.names=F, col.names=T, sep="\t", quote=F)
	write.table(vals, file=paste0(sname,"_vals.txt"), row.names=F, col.names=T, sep="\t", quote=F)

	chr.max <- list()
	for(i in 1:length(chr)) {
		cor.df.chrom <- cor.df[cor.df$chr == chr[i],]
		if (nrow(cor.df.chrom) > 0) {
			chr.max[[i]] <- cor.df.chrom[which.max(cor.df.chrom$score),]
			name <- do.call(rbind,strsplit(clus.df[clus.df$chr == chr[i] & clus.df$group == chr.max[[i]]$group,]$name,"[.]"))
			chr.max[[i]][,"sample"] <- paste0(name[,1], collapse=",")
			chr.max[[i]][,"pcs"]    <- paste0(name[,2], collapse=",")
			chr.max[[i]] <- chr.max[[i]][,-c(ncol(chr.max[[i]])-2)] 
		} else {
			vals.chrom <- vals[vals$chr == as.character(chr[i]),]
			vals.chrom[,"score"] <- as.numeric(vals.chrom$gcc.cor) + as.numeric(vals.chrom$tss.cor)
			sample.df <- list()
			for(j in 1:length(sam)) {
				sample.df[[j]] <- vals.chrom[vals.chrom$sample == sam[j],]	
				sample.df[[j]] <- sample.df[[j]][which.max(sample.df[[j]]$score),]			
			}
			sample.df <- as.data.frame(do.call(rbind, sample.df))
			chr.max[[i]] <- data.frame(
						group=1,
						chr=chr[i],
						gcc.cor=round(mean(sample.df$gcc.cor),4),
						tss.cor=round(mean(sample.df$tss.cor),4),
						len.cor=round(mean(sample.df$len.cor),4),
						sample=paste0(sample.df$sample, collapse=","),
						pcs=paste0(sample.df$pc, collapse=","))	
		}
	}
	chr.max <- (do.call(rbind, chr.max))
	write.table(chr.max, file=paste0(sname,"_chr_pc_selected.txt"), row.names=F, col.names=T, sep="\t", quote=F)
	
	samples <- as.data.frame(do.call(rbind,strsplit(as.character(chr.max$sample),",")))
	pc <- as.data.frame(do.call(rbind,strsplit(as.character(chr.max$pcs),",")))
	for(i in 1:nrow(samples)) {
		for(j in 1:ncol(samples)) {
			cmd <- paste0("cp ",samples[i,j],"_pca/",pc_type,"_pca/",samples[i,j],"_mat/",chr.max$chr[i],".",pc[i,j],".bedGraph ",samples[i,j],"_pca/",pc_type,"_pca/",samples[i,j],"_mat/",chr.max$chr[i],".pc.bedGraph")
			cat ("Running ",cmd,"\n")
			system(cmd, wait=T)
		}
	}	
}

## Select the correct PC and flip 
pcselect <- function(data, genome, pc, diroverwrite, folder=NA) {

	pca_folders <- c("intra_pca","inter_pca")
	
	bed <- read.table(normalizePath(data$bed[1]), h=F, as.is=T)
	resolution <- bed$V3[1] - bed$V2[1]
	chrom_uniq <- unique(bed$V1)

	## Data download and processing 
	if (is.na(folder)) {
    	folder <- paste0(genome,"_",as.integer(resolution),"_goldenpathData")
    	if (!file.exists(folder)) {
    	 	dir.create(folder)
		current_path <- getwd()
        	setwd(folder)
        	datadownload(genome)
        	setwd(current_path)
    	}
    	folder <- normalizePath(folder)
  	} else {
    	folder <- normalizePath(folder)
  	}

  	#current_path <- getwd()
	#setwd(folder)
	#datadownload(genome)
  	#setwd(current_path)

 	if (!file.exists(paste0(folder,"/",genome,".fa"))) {
   		cat ("Unzipping ",paste0(folder,"/",genome,".fa.gz"),"\n")
   		system(paste0("gunzip -c ",folder,"/",genome,".fa.gz > ",folder,"/",genome,".fa"), wait=T)
	}
	if (!file.exists(paste0(folder,"/",genome,".tss.bed"))) {
   		cmd <- paste0("gunzip -c ",folder,"/",genome,".refGene.gtf.gz |awk -v OFS='\\t' '{if($3==\"transcript\"){if($7==\"+\"){print $1,$4,$4+1}else{print $1,$5-1,$5}}}' |grep -v \"alt\" |grep -v \"random\" |sort |uniq |sort -k 1,1 -k2,2n > ",folder,"/",genome,".tss.bed")
   		cat ("Running ",cmd,"\n")
   		system(cmd, wait=T)
	}
	if (!file.exists(paste0(folder,"/",genome,".binned.bed"))) {
  		cmd <- paste0(.findExecutable("bedtools")," makewindows -g ",folder,"/",genome,".chrom.sizes -w ",as.integer(resolution)," > ",folder,"/",genome,".binned.bed")
   		cat ("Running ",cmd,"\n")
   		system(cmd, wait=T)
	}
	if (!file.exists(paste0(folder,"/",genome,".GCpt.bedGraph"))) {
  		cmd <- paste0(.findExecutable("bedtools")," nuc -fi ",folder,"/",genome,".fa -bed ",folder,"/",genome,".binned.bed |grep -v \"#\" |awk -v OFS='\\t' '{print $1,$2,$3,$5}' |grep -v \"alt\" |grep -v \"random\" |sort -k 1,1 -k2,2n > ",folder,"/",genome,".GCpt.bedGraph")
  		cat ("Running ",cmd,"\n")
   		system(cmd, wait=T)
	}
	if (!file.exists(paste0(folder,"/",genome,".GCpt.tss.bedGraph"))) {
		cmd <- paste0(.findExecutable("bedtools")," map -a ",folder,"/",genome,".GCpt.bedGraph -b ",folder,"/",genome,".tss.bed -c 1 -o count -null 0 > ",folder,"/",genome,".GCpt.tss.bedGraph")
    	cat ("Running ",cmd,"\n")
    	system(cmd, wait=T)
  	}

  	gcpt_file <- paste0(folder,"/",genome,".GCpt.tss.bedGraph")
	prefix_master <- unique(data$prefix.master)
	for(n in 1:length(prefix_master)) {

		intra_count <- grep("intra_pca",list.files(paste0(data[data$prefix.master==prefix_master[n],]$prefix,"_pca"), pattern="intra"))
		if (length(intra_count) == nrow(data[data$prefix.master==prefix_master[n],])) {
			pcselectioncore(chrom_uniq, data[data$prefix.master==prefix_master[n],]$prefix, gcpt_file, pc, prefix_master[n], "intra")	
		}

		inter_count <- grep("inter_pca",list.files(paste0(data[data$prefix.master==prefix_master[n],]$prefix,"_pca"), pattern="inter"))
		if (length(inter_count) == nrow(data[data$prefix.master==prefix_master[n],])) {
			pcselectioncore(chrom_uniq, data[data$prefix.master==prefix_master[n],]$prefix, gcpt_file, pc, prefix_master[n], "inter")	
		}
	}
}

### END ###

########################### Get differential results from PC data ###########################

options(width=400)
#Mahalabobis weight
mdweight <- function(m1) {
  d1 <- list()
  w1 <- list()
  for(i in 1:length(m1)) {
    d1[[i]] <- sum(abs(m1[i] - m1))
    w1[[i]] <- (length(m1) * d1[[i]])
  }
  d1 <- unlist(d1)
  w1 <- unlist(w1)
  w1 <- w1/ifelse(sum(d1)==0, 1, sum(d1))
  return(list(weight=w1))
}

#Calculate centers
calcen <- function(df, class, rzscore, szscore, chrom) {

	df_dist <- as.data.frame(t(apply(df, 1, function(x){sqrt(colSums(as.matrix(dist(as.numeric(x))))/(length(as.numeric(x))-1))})))
	colnames(df_dist) <- paste0("sqrtdist.",colnames(df))
	
	#Zscore estimate
	if (class == "rep") {
		value <- as.numeric(as.matrix(df_dist))
		col_mean  <- mean(value)
		col_sd    <- sd(value)
		df_zsc 	  <- list()
		for(n in 1:ncol(df_dist)) {
                	df_zsc[[n]] <- (df_dist[,n] - col_mean)/col_sd
		}
		df_zsc <- do.call(cbind,df_zsc)
		df_pvl <- round(pnorm(df_zsc, mean=rzscore, lower.tail=T),5)
		df_pvl_max <- apply(df_pvl, 1, max)
		df_cen <- round((df * (1 - df_pvl_max)),5)
		return(df_cen)

	} else if (class == "sample") { 
		df_zsc <- list()
		for(n in 1:ncol(df)) {
			value <- as.numeric(as.matrix(df_dist[,n]))
			col_mean  <- mean(value)
	        col_sd    <- sd(value)
			df_zsc[[n]] <- (value - col_mean)/col_sd
		}
		df_zsc <- do.call(cbind,df_zsc)
		colnames(df_zsc) <- paste0("zsc.",colnames(df))
		df_pvl <- round(pnorm(as.matrix(df_zsc), mean=szscore, lower.tail=T),5)
		df_pvl_max <- apply(df_pvl, 1, max)
		df_cen <- round((df * (1 - df_pvl_max)),5)
		colnames(df_cen) <- paste0("samcen.",colnames(df))
		return(df_cen)
	}
}

calcGlosh <- function(df, fdr.thr, glosh_k=5) {

	cat("Calculating Glosh score\n")
	rownames(df) <- paste0(df$chr,"_",df$start)
	df_sig <- df[df$padj < fdr.thr,]
	chrom  <- as.character(unique(df_sig$chr))
	gloshscore <- list()
	for(i in 1:length(chrom)) {
		gloshscore[[i]] <- dbscan::glosh(as.matrix(dist(df_sig[df_sig$chr==as.character(chrom[i]),]$start)/1e6), k=glosh_k)
	}
	gloshscore <- unlist(gloshscore)
	df_sig[,"glosh"] <- gloshscore
	df[,"glosh"] <- 1
	df[rownames(df_sig),"glosh"] <- df_sig$glosh
	df$glosh <- 1 - df$glosh
	return(df)
}


calchclust <- function(df, fdr.thr, hclust_dist=distclust) {

	cat("Clustering by distance and removing lone peaks\n")
	rownames(df) <- paste0(df$chr,"_",df$start)
	df_sig <- df[df$padj < fdr.thr,]
	chrom  <- as.character(unique(df_sig$chr))
	dist_clust <- list()
	for(i in 1:length(chrom)) {
		cl <- cutree(hclust(dist(df_sig[df_sig$chr==as.character(chrom[i]),]$start)/1e6), h=hclust_dist/1e6)
		dist_clust[[i]] <- as.integer(table(cl)[cl])
	}
	dist_clust <- unlist(dist_clust)
	df_sig[,"dist_clust"] <- dist_clust
	df[,"dist_clust"] <- 1
	df[rownames(df_sig),"dist_clust"] <- df_sig$dist_clust
	return(df)
}


pcanalyze <- function(data, diroverwrite, diffolder, rzscore, szscore, fdr.thr, distclust, numberclust) {

	#Reformatting the data
	if (!dir.exists("DifferentialResult") | diroverwrite == TRUE) {
		if (!dir.exists("DifferentialResult")) {
			dir.create("DifferentialResult")
		}
		if (!dir.exists(paste0("DifferentialResult/",diffolder))) {
			dir.create(paste0("DifferentialResult/",diffolder))
		}

		diffolder    <- paste0("DifferentialResult/",diffolder)
		intra_bed_df <- list()
		inter_bed_df <- list()
		intra_counter<- 0
		inter_counter<- 0
		sample_count <- 0
		prefix_master <- unique(data$prefix.master)
		for(i in 1:length(prefix_master)) {
			
			data_rep <- data[data$prefix.master == prefix_master[i],]
			sample_count <- sample_count + nrow(data_rep)
			if (!dir.exists(paste0(diffolder,"/",prefix_master[i],"_data"))) {
				dir.create(paste0(diffolder,"/",prefix_master[i],"_data"))
			}

			for(j in 1:nrow(data_rep)) {
				if (dir.exists(paste0(data_rep$prefix[j],"_pca/intra_pca/",data_rep$prefix[j],"_mat"))) {
					chrom <- list.files(paste0(data_rep$prefix[j],"_pca/intra_pca/",data_rep$prefix[j],"_mat"), pattern=".pc.bedGraph", full.names=T)
					chrom <- gsub(".pc.bedGraph","",chrom)
					for(k in 1:length(chrom)) {
						cmd <- paste0("cp ",chrom[k],".pc.bedGraph ",diffolder,"/",prefix_master[i],"_data/",data_rep$prefix[j],"_intra_",basename(chrom[k]),".pc.bedGraph")
						system(cmd, wait=T)
					}
				}
				if (dir.exists(paste0(data_rep$prefix[j],"_pca/inter_pca/",data_rep$prefix[j],"_mat"))) {
					chrom <- list.files(paste0(data_rep$prefix[j],"_pca/inter_pca/",data_rep$prefix[j],"_mat"), pattern=".pc.bedGraph", full.names=T)
					chrom <- gsub(".pc.bedGraph","",chrom)
					for(k in 1:length(chrom)) {
						cmd <- paste0("cp ",chrom[k],".pc.bedGraph ",diffolder,"/",prefix_master[i],"_data/",data_rep$prefix[j],"_inter_",basename(chrom[k]),".pc.bedGraph")
						system(cmd, wait=T)
					}
				}
			}

			for(k in 1:length(chrom)) {
				chrom_rep <- list.files(paste0(diffolder,"/",prefix_master[i],"_data"), pattern=paste0("_intra_",basename(chrom[k]),".pc.bedGraph"), full.names=T)
				if (length(chrom_rep) > 0) {
					pcmat   <- list()
					rnames  <- list()
					namevec <- c()
					for(u in 1:length(chrom_rep)) {
						df <- read.table(chrom_rep[u], h=F, as.is=T)
						colnames(df)<- c("chr","start","end","PC")
						rownames(df)<- paste0(df$chr,"_",df$start,"_",df$end)
						namevec[u]  <- strsplit(gsub("_intra"," ",basename(chrom_rep[u]))," ")[[1]][1]
						rnames[[u]] <- df[,1:3]
						intra_counter <- intra_counter + 1
						intra_bed_df[[intra_counter]] <- df[,1:3]
						pcmat[[u]] <- df
					}
					rnames <- do.call(rbind, rnames)
					rnames <- as.data.frame(data.table::setDT(rnames)[,list(Count=.N),names(rnames)])
					rnames <- rnames[rnames$Count == length(chrom_rep),1:3]
					bed    <- rnames
					rnames <- paste0(rnames$chr,"_",rnames$start,"_",rnames$end)
					for(u in 1:length(pcmat)) {
						pcmat[[u]] <- pcmat[[u]][as.character(rnames),4]
					}
					pcmat    <- do.call(cbind, pcmat)
					pcmat_qn <- limma::normalizeQuantiles(pcmat, ties=T)
					pcmat_qn <- pcmat_qn/apply(abs(pcmat_qn), 2, mean)
					bed      <- cbind(bed,pcmat)
					bed_qn   <- cbind(bed[,1:3],pcmat_qn)
					colnames(bed)    <- c("chr","start","end",namevec)
					colnames(bed_qn) <- c("chr","start","end",namevec)
					write.table(bed, file=paste0(diffolder,"/",prefix_master[i],"_data/intra_",basename(chrom[k]),"_combined.pcOri.bedGraph"), row.names=F, col.names=T, sep="\t", quote=F)
					write.table(bed_qn, file=paste0(diffolder,"/",prefix_master[i],"_data/intra_",basename(chrom[k]),"_combined.pcQnm.bedGraph"), row.names=F, col.names=T, sep="\t", quote=F)
					cat("Wrote ",paste0(diffolder,"/",prefix_master[i],"_data/intra_",basename(chrom[k]),"_combined.pcOri.bedGraph and pcQnm.bedGraph files\n"))
				}

				chrom_rep <- list.files(paste0(diffolder,"/",prefix_master[i],"_data"), pattern=paste0("_inter_",basename(chrom[k]),".pc.bedGraph"), full.names=T)
				if (length(chrom_rep) > 0) {
					pcmat   <- list()
					rnames  <- list()
					namevec <- c()
					for(u in 1:length(chrom_rep)) {
						df <- read.table(chrom_rep[u], h=F, as.is=T)
						colnames(df)<- c("chr","start","end","PC")
						rownames(df)<- paste0(df$chr,"_",df$start,"_",df$end)
						namevec[u]  <- strsplit(gsub("_inter"," ",basename(chrom_rep[u]))," ")[[1]][1]
						rnames[[u]] <- df[,1:3]
						inter_counter <- inter_counter + 1
						inter_bed_df[[inter_counter]] <- df[,1:3]
						pcmat[[u]] <- df
					}
					rnames <- do.call(rbind, rnames)
					rnames <- as.data.frame(data.table::setDT(rnames)[,list(Count=.N),names(rnames)])
					rnames <- rnames[rnames$Count == length(chrom_rep),1:3]
					bed    <- rnames
					rnames <- paste0(rnames$chr,"_",rnames$start,"_",rnames$end)
					for(u in 1:length(pcmat)) {
						pcmat[[u]] <- pcmat[[u]][as.character(rnames),4]
					}
					pcmat 	 <- do.call(cbind, pcmat)
					pcmat_qn <- limma::normalizeQuantiles(pcmat, ties=T)
					pcmat_qn <- pcmat_qn/apply(abs(pcmat_qn), 2, mean)
					bed      <- cbind(bed,pcmat)
					bed_qn   <- cbind(bed[,1:3],pcmat_qn)
					colnames(bed)    <- c("chr","start","end",namevec)
					colnames(bed_qn) <- c("chr","start","end",namevec)
					write.table(bed, file=paste0(diffolder,"/",prefix_master[i],"_data/inter_",basename(chrom[k]),"_combined.pcOri.bedGraph"), row.names=F, col.names=T, sep="\t", quote=F)
					write.table(bed_qn, file=paste0(diffolder,"/",prefix_master[i],"_data/inter_",basename(chrom[k]),"_combined.pcQnm.bedGraph"), row.names=F, col.names=T, sep="\t", quote=F)
					cat("Wrote ",paste0(diffolder,"/",prefix_master[i],"_data/inter_",basename(chrom[k]),"_combined.pcOri.bedGraph and pcQnm.bedGraph files\n"))
				}
			}
		}

		if (!dir.exists(paste0(diffolder,"/pcOri"))) {
			dir.create(paste0(diffolder,"/pcOri"))
		}
		if (!dir.exists(paste0(diffolder,"/pcQnm"))) {
			dir.create(paste0(diffolder,"/pcQnm"))
		}

		if (intra_counter > 0) {
			intra_bed_df <- as.data.frame(do.call(rbind, intra_bed_df))
			colnames(intra_bed_df) <- c("chr","start","end")
			intra_bed_df <- as.data.frame(data.table::setDT(intra_bed_df)[,list(Count=.N),names(intra_bed_df)])
		}

		if (inter_counter > 0) {
			inter_bed_df <- as.data.frame(do.call(rbind, inter_bed_df))
			colnames(inter_bed_df) <- c("chr","start","end")
			inter_bed_df <- as.data.frame(data.table::setDT(inter_bed_df)[,list(Count=.N),names(inter_bed_df)])
		}

		chrom <- basename(chrom)
		for(k in 1:length(chrom)) {

			if (intra_counter > 0) {
				chrom_sample <- list()
				namevec <- c()
				intra_bed_chr <- intra_bed_df[intra_bed_df$chr == as.character(chrom[k]) & intra_bed_df$Count == sample_count,]
				intra_bed_vec <- paste0(intra_bed_chr$chr,"_",intra_bed_chr$start,"_",intra_bed_chr$end) 
				for(i in 1:length(prefix_master)) {
					if (file.exists(paste0(diffolder,"/",prefix_master[i],"_data/intra_",chrom[k],"_combined.pcOri.bedGraph"))) {
						df <- read.table(paste0(diffolder,"/",prefix_master[i],"_data/intra_",chrom[k],"_combined.pcOri.bedGraph"), h=T, as.is=T)
						rownames(df) <- paste0(df$chr,"_",df$start,"_",df$end)
						if (i == 1) {
							chrom_sample[[i]] <- df[as.character(intra_bed_vec),]
							namevec <- c(namevec, colnames(df))
						} else {
							chrom_sample[[i]] <- df[as.character(intra_bed_vec),c(4:ncol(df))]
							namevec <- c(namevec, colnames(df)[4:ncol(df)])
						}
					}
				}
				if (length(chrom_sample) > 0) {
					chrom_sample    <- as.data.frame(do.call(cbind, chrom_sample))
					chrom_sample_qn <- limma::normalizeQuantiles(chrom_sample[,4:ncol(chrom_sample)], ties=T)
					chrom_sample_qn <- chrom_sample_qn/apply(abs(chrom_sample_qn), 2, mean)
					chrom_sample_qn <- cbind(chrom_sample[,1:3],chrom_sample_qn)
					colnames(chrom_sample)    <- namevec
					colnames(chrom_sample_qn) <- namevec
					chrom_sample[,4:ncol(chrom_sample)] <- round(chrom_sample[,4:ncol(chrom_sample)], 5)
					chrom_sample_qn[,4:ncol(chrom_sample_qn)] <- round(chrom_sample_qn[,4:ncol(chrom_sample_qn)], 5)
					write.table(chrom_sample, file=paste0(diffolder,"/pcOri/intra_sample_",chrom[k],"_combined.pcOri.bedGraph"), row.names=F, sep="\t", quote=F)
					write.table(chrom_sample_qn, file=paste0(diffolder,"/pcQnm/intra_sample_",chrom[k],"_combined.pcQnm.bedGraph"), row.names=F, sep="\t", quote=F)
					cat (paste0("Wrote intra_sample_",chrom[k],"_combined.pcOri.bedGraph & _combined.pcQnm.bedGraph files under ",diffolder,"/pcOri & ",diffolder,"/pcQnm folders"),"\n")
				}
			}

			if (inter_counter > 0) {
				chrom_sample <- list()
				namevec <- c()
				inter_bed_chr <- inter_bed_df[inter_bed_df$chr == as.character(chrom[k]) & inter_bed_df$Count == sample_count,]
				inter_bed_vec <- paste0(inter_bed_chr$chr,"_",inter_bed_chr$start,"_",inter_bed_chr$end) 
				for(i in 1:length(prefix_master)) {
					if (file.exists(paste0(diffolder,"/",prefix_master[i],"_data/inter_",chrom[k],"_combined.pcOri.bedGraph"))) {
						df <- read.table(paste0(diffolder,"/",prefix_master[i],"_data/inter_",chrom[k],"_combined.pcOri.bedGraph"), h=T, as.is=T)
						rownames(df) <- paste0(df$chr,"_",df$start,"_",df$end)
						if (i == 1) {
							chrom_sample[[i]] <- df[as.character(inter_bed_vec),]
							namevec <- c(namevec, colnames(df))
						} else {
							chrom_sample[[i]] <- df[as.character(inter_bed_vec),c(4:ncol(df))]
							namevec <- c(namevec, colnames(df)[4:ncol(df)])
						}
					}
				}
				if (length(chrom_sample) > 0) {
					chrom_sample    <- do.call(cbind, chrom_sample)
					chrom_sample_qn <- limma::normalizeQuantiles(chrom_sample[,4:ncol(chrom_sample)], ties=T)
					chrom_sample_qn <- chrom_sample_qn/apply(abs(chrom_sample_qn), 2, mean)
					chrom_sample_qn <- cbind(chrom_sample[,1:3],chrom_sample_qn)
					colnames(chrom_sample)    <- namevec
					colnames(chrom_sample_qn) <- namevec
					chrom_sample[,4:ncol(chrom_sample)] <- round(chrom_sample[,4:ncol(chrom_sample)], 5)
					chrom_sample_qn[,4:ncol(chrom_sample_qn)] <- round(chrom_sample_qn[,4:ncol(chrom_sample_qn)], 5)
					write.table(chrom_sample, file=paste0(diffolder,"/pcOri/inter_sample_",chrom[k],"_combined.pcOri.bedGraph"), row.names=F, sep="\t", quote=F)
					write.table(chrom_sample_qn, file=paste0(diffolder,"/pcQnm/inter_sample_",chrom[k],"_combined.pcQnm.bedGraph"), row.names=F, sep="\t", quote=F)
					cat (paste0("Wrote inter_sample_",chrom[k],"_combined.pcOri.bedGraph & _combined.pcQnm.bedGraph files under ",diffolder,"/pcOri & ",diffolder,"/pcQnm folders"),"\n")
				}
			}
		}

		#Intra significance calculation
		intra_bdg <- list.files(paste0(diffolder,"/pcQnm"),pattern="intra",full.names=T)
		intra_bdg_pcori <- list.files(paste0(diffolder,"/pcOri"),pattern="intra",full.names=T)
		intra_res <- list()
		grp_pcOri <- list()
		if (length(intra_bdg) > 0) {
			pval_wt <- list()
			for(j in 1:length(intra_bdg)) {
				intra_cen <- list()
				intra_grp <- list()
				intra_grp_pcori <- list()
				rep_count <- 0
				df_intra  <- read.table(intra_bdg[j], h=T, as.is=T)
				df_intra_pcori <- read.table(intra_bdg_pcori[j], h=T, as.is=T) 
				for(i in 1:length(prefix_master)) {
					data_rep <- data[data$prefix.master == prefix_master[i],]
					if (nrow(data_rep) > 1) {
						intra_cen[[i]] <- calcen(df_intra[,data_rep$prefix], class="rep", rzscore, szscore, df_intra$chr[1])
						intra_grp[[i]] <- apply(df_intra[,data_rep$prefix], 1, mean)
						intra_grp_pcori[[i]] <- apply(df_intra_pcori[,data_rep$prefix], 1, mean) 
						rep_count <- rep_count + 1
					} else {
						intra_cen[[i]] <- rep(mean(df_intra[,data_rep$prefix]), nrow(df_intra))
						intra_grp[[i]] <- df_intra[,data_rep$prefix]
						intra_grp_pcori[[i]] <- df_intra_pcori[,data_rep$prefix]
					}
				}
				intra_cen <- do.call(cbind, intra_cen)
				intra_grp <- do.call(cbind, intra_grp)
				grp_pcOri[[j]] <- as.data.frame(do.call(cbind, intra_grp_pcori))
				colnames(intra_cen) <- paste0("repcen.",colnames(df_intra)[-c(1:3)])
				colnames(intra_grp) <- prefix_master
				colnames(grp_pcOri[[j]]) <- prefix_master
				rownames(grp_pcOri[[j]]) <- paste0(df_intra_pcori$chr,"_",df_intra_pcori$start) 

				if (rep_count > 0) {
					replicate_maha <- list()
					maha_temp <- list()
					for(k in 1:nrow(df_intra)) {
						replicate_maha[[k]] <- mahalanobis(as.matrix(df_intra[k,c(4:ncol(df_intra))]), as.matrix(intra_cen[k,]), diag(diag(cov(df_intra[,c(4:ncol(df_intra))]))))
					}
					replicate_maha <- unlist(replicate_maha)
					replicate_wt   <- -pchisq(replicate_maha, df=ncol(df_intra)-4, lower.tail=T, log.p=T)
				}

				intra_grp_cen <- calcen(intra_grp, class="sample", rzscore, szscore, df_intra$chr[1])
				max_min_dif <- abs(apply(intra_grp, 1, max) - apply(intra_grp, 1, min))
				sam_mcd.cov <- robust::covRob(intra_grp)$cov
				sample_maha <- list()
				for(k in 1:nrow(intra_grp)) {
					sample_maha[[k]] <- mahalanobis(intra_grp[k,], as.matrix(intra_grp_cen[k,]), (sam_mcd.cov * (1/max_min_dif[k]^2)))
				}
				
				sample_maha <- unlist(sample_maha) 
				pval <- pchisq(sample_maha, df=ncol(intra_grp)-1, lower.tail=F)

				if (rep_count > 0) {
					intra_res[[j]]  <- data.frame(df_intra,intra_grp,replicate_wt,sample_maha,pval)
				} else {
					intra_res[[j]]  <- data.frame(df_intra,intra_grp,sample_maha,pval)
				}
				
				if (!dir.exists(paste0(diffolder,"/fdr_result"))) {
					dir.create(paste0(diffolder,"/fdr_result"))
				}
				write.table(intra_res[[j]], file=paste0(diffolder,"/fdr_result/differential.",basename(intra_bdg[j])), row.names=F, sep="\t", quote=F)
				cat ("Wrote ",paste0(diffolder,"/fdr_result/differential.",basename(intra_bdg[j]))," file\n")
			}
			if(length(intra_res) > 0) { 
				intra_res <- as.data.frame(do.call(rbind,intra_res))
				intra_res <- na.omit(intra_res)
				if (rep_count > 0) {
					ihw_val <- IHW::ihw(pval ~ replicate_wt,  data = intra_res, alpha = 0.2)
					intra_res[,"padj"] <- IHW::adj_pvalues(ihw_val)	
					intra_res <- intra_res[order(intra_res$chr, intra_res$start),]
				} else {
					intra_res[,"padj"] <- p.adjust(intra_res$pval, "BH")
				}

				if (distclust > -1) {
					#intra_res <- calcGlosh(intra_res, fdr.thr, glosh_k=5)
					intra_res <- calchclust(intra_res, fdr.thr, hclust_dist=distclust)
				} else if (distclust == -1) {
					#intra_res[,"glosh"] <- 1
					intra_res[,"dist_clust"] <- 2
					rownames(intra_res) <- paste0(intra_res$chr,"_",intra_res$start)
				}

				write.table(intra_res, file=paste0(diffolder,"/fdr_result/differential.intra_sample_combined.pcQnm.bedGraph"), row.names=F, sep="\t", quote=F)
				cat ("Wrote ",paste0(diffolder,"/fdr_result/differential.intra_sample_combined.pcQnm.bedGraph")," file\n")

				if (fdr.thr > 0) {
					#write.table(intra_res[intra_res$padj < fdr.thr & intra_res$glosh >= gscore, ], file=paste0(diffolder,"/fdr_result/differential.intra_sample_combined.Filtered.pcQnm.bedGraph"), row.names=F, sep="\t", quote=F)
					write.table(intra_res[intra_res$padj < fdr.thr & intra_res$dist_clust > numberclust, ], file=paste0(diffolder,"/fdr_result/differential.intra_sample_combined.Filtered.pcQnm.bedGraph"), row.names=F, sep="\t", quote=F)
				}

				intra_res <- intra_res[,-c(which(colnames(intra_res) %in% data$prefix))]
				write.table(intra_res, file=paste0(diffolder,"/fdr_result/differential.intra_sample_group.pcQnm.bedGraph"), row.names=F, sep="\t", quote=F)
				cat ("Wrote ",paste0(diffolder,"/fdr_result/differential.intra_sample_group.pcQnm.bedGraph")," file\n")

				if (fdr.thr > 0) {
					#write.table(intra_res[intra_res$padj < fdr.thr & intra_res$glosh >= gscore, ], file=paste0(diffolder,"/fdr_result/differential.intra_sample_group.Filtered.pcQnm.bedGraph"), row.names=F, sep="\t", quote=F)
					write.table(intra_res[intra_res$padj < fdr.thr & intra_res$dist_clust > numberclust, ], file=paste0(diffolder,"/fdr_result/differential.intra_sample_group.Filtered.pcQnm.bedGraph"), row.names=F, sep="\t", quote=F)
				}

				grp_pcOri <- as.data.frame(do.call(rbind, grp_pcOri))
				intra_res[,as.character(unique(data$prefix.master))]  <- grp_pcOri[rownames(intra_res),as.character(unique(data$prefix.master))]
				write.table(intra_res, file=paste0(diffolder,"/fdr_result/differential.intra_sample_group.pcOri.bedGraph"), row.names=F, sep="\t", quote=F)
				cat ("Wrote ",paste0(diffolder,"/fdr_result/differential.intra_sample_group.pcOri.bedGraph")," file\n")

				if (fdr.thr > 0) {
					#write.table(intra_res[intra_res$padj < fdr.thr & intra_res$glosh >= gscore, ], file=paste0(diffolder,"/fdr_result/differential.intra_sample_group.Filtered.pcOri.bedGraph"), row.names=F, sep="\t", quote=F)
					write.table(intra_res[intra_res$padj < fdr.thr & intra_res$dist_clust > numberclust, ], file=paste0(diffolder,"/fdr_result/differential.intra_sample_group.Filtered.pcOri.bedGraph"), row.names=F, sep="\t", quote=F)
				}
			}
		}

		#Inter significance calculation
		inter_bdg <- list.files(paste0(diffolder,"/pcQnm"),pattern="inter",full.names=T)
		inter_bdg_pcori <- list.files(paste0(diffolder,"/pcOri"),pattern="inter",full.names=T)
		grp_pcOri <- list()
		inter_res <- list()
		if (length(inter_bdg) > 0) {
			pval_wt <- list()
			for(j in 1:length(inter_bdg)) {
				inter_cen <- list()
				inter_grp <- list()
				inter_grp_pcori <- list()
				rep_count <- 0
				df_inter  <- read.table(inter_bdg[j], h=T, as.is=T)
				df_inter_pcori <- read.table(inter_bdg_pcori[j], h=T, as.is=T) 
				for(i in 1:length(prefix_master)) {
					data_rep <- data[data$prefix.master == prefix_master[i],]
					if (nrow(data_rep) > 1) {
						inter_cen[[i]] <- calcen(df_inter[,data_rep$prefix], class="rep", rzscore, szscore, df_inter$chr[1])
						inter_grp[[i]] <- apply(df_inter[,data_rep$prefix], 1, mean)
						inter_grp_pcori[[i]] <- apply(df_inter_pcori[,data_rep$prefix], 1, mean) 
						rep_count <- rep_count + 1
					} else {
						inter_cen[[i]] <- rep(mean(df_inter[,data_rep$prefix]), nrow(df_inter))
						inter_grp[[i]] <- df_inter[,data_rep$prefix]
						inter_grp_pcori[[i]] <- df_inter_pcori[,data_rep$prefix]
					}
				}
				inter_cen <- do.call(cbind, inter_cen)
				inter_grp <- do.call(cbind, inter_grp)
				grp_pcOri[[j]] <- as.data.frame(do.call(cbind, inter_grp_pcori))
				colnames(inter_cen) <- paste0("repcen.",colnames(df_inter)[-c(1:3)])
				colnames(inter_grp) <- prefix_master
				colnames(grp_pcOri[[j]]) <- prefix_master
				rownames(grp_pcOri[[j]]) <- paste0(df_inter_pcori$chr,"_",df_inter_pcori$start) 

				if (rep_count > 0) {
					replicate_maha <- list()
					maha_temp <- list()
					for(k in 1:nrow(df_inter)) {
						replicate_maha[[k]] <- mahalanobis(as.matrix(df_inter[k,c(4:ncol(df_inter))]), as.matrix(inter_cen[k,]), diag(diag(cov(df_inter[,c(4:ncol(df_inter))]))))
					}
					replicate_maha <- unlist(replicate_maha)
					replicate_wt   <- -pchisq(replicate_maha, df=ncol(df_inter)-4, lower.tail=T, log.p=T)
				}
				
				inter_grp_cen <- calcen(inter_grp, class="sample", rzscore, szscore, df_inter$chr[1])
				max_min_dif <- abs(apply(inter_grp, 1, max) - apply(inter_grp, 1, min))
				sam_mcd.cov <- robust::covRob(inter_grp)$cov
				sample_maha <- list()
				for(k in 1:nrow(inter_grp)) {
					sample_maha[[k]] <- mahalanobis(inter_grp[k,], as.matrix(inter_grp_cen[k,]), (sam_mcd.cov * (1/max_min_dif[k]^2)))
				}
				
				sample_maha <- unlist(sample_maha) 
				pval <- pchisq(sample_maha, df=ncol(inter_grp)-1, lower.tail=F)

				if (rep_count > 0) {
					inter_res[[j]]  <- data.frame(df_inter,inter_grp,replicate_wt,sample_maha,pval)
				} else {
					inter_res[[j]]  <- data.frame(df_inter,inter_grp,sample_maha,pval)
				}
				
				if (!dir.exists(paste0(diffolder,"/fdr_result"))) {
					dir.create(paste0(diffolder,"/fdr_result"))
				}
				write.table(inter_res[[j]], file=paste0(diffolder,"/fdr_result/differential.",basename(inter_bdg[j])), row.names=F, sep="\t", quote=F)
				cat ("Wrote ",paste0(diffolder,"/fdr_result/differential.",basename(inter_bdg[j]))," file\n")
			}
			if(length(inter_res) > 0) { 
				inter_res <- as.data.frame(do.call(rbind,inter_res))
				inter_res <- na.omit(inter_res)
				if (rep_count > 0) {
					ihw_val <- IHW::ihw(pval ~ replicate_wt,  data = inter_res, alpha = 0.2)
					inter_res[,"padj"] <- IHW::adj_pvalues(ihw_val)	
					inter_res <- inter_res[order(inter_res$chr, inter_res$start),]
				} else {
					inter_res[,"padj"] <- p.adjust(inter_res$pval, "BH")
				}

				if (distclust > -1) {
					#inter_res <- calcGlosh(inter_res, fdr.thr, glosh_k=5)
					inter_res <- calchclust(inter_res, fdr.thr, hclust_dist=distclust)
				} else if (distclust == -1) {
					#inter_res[,"glosh"] <- 1
					inter_res[,"dist_clust"] <- 2 
					rownames(inter_res) <- paste0(inter_res$chr,"_",inter_res$start)
				}

				write.table(inter_res, file=paste0(diffolder,"/fdr_result/differential.inter_sample_combined.pcQnm.bedGraph"), row.names=F, sep="\t", quote=F)
				cat ("Wrote ",paste0(diffolder,"/fdr_result/differential.inter_sample_combined.pcQnm.bedGraph")," file\n")

				if (fdr.thr > 0) {
					#write.table(inter_res[inter_res$padj < fdr.thr & inter_res$glosh >= gscore, ], file=paste0(diffolder,"/fdr_result/differential.inter_sample_combined.Filtered.pcQnm.bedGraph"), row.names=F, sep="\t", quote=F)
					write.table(inter_res[inter_res$padj < fdr.thr & inter_res$dist_clust > numberclust, ], file=paste0(diffolder,"/fdr_result/differential.inter_sample_combined.Filtered.pcQnm.bedGraph"), row.names=F, sep="\t", quote=F)
				}

				inter_res <- inter_res[,-c(which(colnames(inter_res) %in% data$prefix))]
				write.table(inter_res, file=paste0(diffolder,"/fdr_result/differential.inter_sample_group.pcQnm.bedGraph"), row.names=F, sep="\t", quote=F)
				cat ("Wrote ",paste0(diffolder,"/fdr_result/differential.inter_sample_group.pcQnm.bedGraph")," file\n")

				if (fdr.thr > 0) {
					#write.table(inter_res[inter_res$padj < fdr.thr & inter_res$glosh >= gscore, ], file=paste0(diffolder,"/fdr_result/differential.inter_sample_group.Filtered.pcQnm.bedGraph"), row.names=F, sep="\t", quote=F)
					write.table(inter_res[inter_res$padj < fdr.thr & inter_res$dist_clust > numberclust, ], file=paste0(diffolder,"/fdr_result/differential.inter_sample_group.Filtered.pcQnm.bedGraph"), row.names=F, sep="\t", quote=F)
				}

				grp_pcOri <- as.data.frame(do.call(rbind, grp_pcOri))
				inter_res[,as.character(unique(data$prefix.master))]  <- grp_pcOri[rownames(inter_res),as.character(unique(data$prefix.master))]
				write.table(inter_res, file=paste0(diffolder,"/fdr_result/differential.inter_sample_group.pcOri.bedGraph"), row.names=F, sep="\t", quote=F)
				cat ("Wrote ",paste0(diffolder,"/fdr_result/differential.inter_sample_group.pcOri.bedGraph")," file\n")

				if (fdr.thr > 0) {
					#write.table(inter_res[inter_res$padj < fdr.thr & inter_res$glosh >= gscore, ], file=paste0(diffolder,"/fdr_result/differential.inter_sample_group.Filtered.pcOri.bedGraph"), row.names=F, sep="\t", quote=F)
					write.table(inter_res[inter_res$padj < fdr.thr & inter_res$dist_clust > numberclust, ], file=paste0(diffolder,"/fdr_result/differential.inter_sample_group.Filtered.pcOri.bedGraph"), row.names=F, sep="\t", quote=F)
				}
			}
		}
	}
}

### END ###

########################### Sub-compartments finding ###########################

subcompartment_level <- function(vec) {
	
	A_comp <- vec[vec > 0]
	B_comp <- vec[vec < 0]

	A_list <- list()
	B_list <- list()

	B_list[[1]] <- data.frame(score1=B_comp[1], score2=B_comp[2], mean_score=mean(c(B_comp[1],B_comp[2])), label=paste0("B",(length(B_comp)-1)))
	i <- 2
	while (i < length(B_comp)) {
		B_list[[i]] <- data.frame(score1=B_comp[i], score2=B_comp[(i+1)], mean_score=mean(c(B_comp[i],B_comp[(i+1)])), label=paste0("B",(length(B_comp)-i)))
		i <- i + 1
	}
	B_list[[i]] <- data.frame(score1=B_comp[i], score2=0, mean_score=mean(c(B_comp[i],0)), label=paste0("B0"))

	A_list[[1]] <- data.frame(score1=0, score2=A_comp[1], mean_score=mean(c(0,A_comp[1])), label=paste0("A0"))
	i <- 2
	while (i <= length(A_comp)) {
		A_list[[i]] <- data.frame(score1=A_comp[(i-1)], score2=A_comp[i], mean_score=mean(c(A_comp[(i-1)],A_comp[i])), label=paste0("A",(i-1)))
		i <- i + 1
	}

	B_list <- do.call(rbind, B_list)
	A_list <- do.call(rbind, A_list)
	df <- as.data.frame(rbind(B_list,A_list))
	return(df)
}

hmmsegment <- function(compartment_file, prefix_master, subnum) {

	chrom  <- unique(compartment_file$chr)
	hmm_comp <- list()
	for(i in 1:length(chrom)) {
		cat ("Running hmm segmentation on ",chrom[i],"\n")
		chrom_file  <- compartment_file[compartment_file$chr==as.character(chrom[i]),]
		state_score <- list()
		labels <- subcompartment_level(as.vector((quantile(c(as.matrix(chrom_file[,c(prefix_master)])),probs=c(0:subnum)/subnum))))
		
		for(j in 1:length(prefix_master)) {
			cat ("Runnig for ",prefix_master[j],"\n")
			df <- chrom_file[,c(1:3,grep(prefix_master[j],colnames(chrom_file)))]
			colnames(df) <- c("chr","start","end","sample")
			mod_f  <- depmixS4::depmix(sample ~ 1, data=df, nstates=subnum, family=gaussian())
			mod_b  <- depmixS4::depmix(sample ~ 1, data=df[order(-df$start),], nstates=subnum, family=gaussian())
			cat ("Running forward segmentation\n")
			fmod_f <- depmixS4::fit(mod_f)
			cat ("Running backward segmentation\n")
			fmod_b <- depmixS4::fit(mod_b)
			df[,"state_f"] <- fmod_f@posterior$state
			df[,"state_b"] <- rev(fmod_b@posterior$state)
			print (head(df))
			state_f <- aggregate(sample ~ state_f, mean, data=df)
			state_b <- aggregate(sample ~ state_b, mean, data=df)
			state_f <- state_f[order(state_f$sample),]
			state_b <- state_b[order(state_b$sample),]
			state_f[,"state"] <- c(1:subnum)
			state_b[,"state"] <- c(1:subnum)
			print (state_f)
			print (state_b)
			df[,"state"] <- 0
			for(k in 1:nrow(df)) {
				old_s_f <- df$state_f[k]
				old_s_b <- df$state_b[k]
				new_s_f <- state_f[state_f$state_f == old_s_f, "state"]
				new_s_b <- state_b[state_b$state_b == old_s_b, "state"]
				s <- ceiling(mean(c(new_s_f,new_s_b)))
				df$state[k] <- s
			}
			print (head(df))
			state_score[[j]] <- aggregate(sample ~ state, mean, data=df)
			chrom_file[,"state"] <- df$state
			colnames(chrom_file)[ncol(chrom_file)] <- paste0(prefix_master[j],".state")
			chrom_file[,"score"] <- df$state
			colnames(chrom_file)[ncol(chrom_file)] <- paste0(prefix_master[j],".score")
		}
		
		chrom_comp <- list()
		for(j in 1:length(state_score)) {
			state_score[[j]][,"label"] <- "XX"
			state_score[[j]][,"score"] <- 0
			for(k in 1:nrow(labels)) {
				for(l in 1:nrow(state_score[[j]])) {
					if (state_score[[j]]$sample[l] >= labels$score1[k] & state_score[[j]]$sample[l] < labels$score2[k]) {
						state_score[[j]][l,"label"] <- as.character(labels$label[k])
						state_score[[j]][l,"score"] <- round(labels$mean_score[k],3)
					}
				}
			}
			pos_label <- grep(paste0(prefix_master[j],".state"), colnames(chrom_file))
			pos_score <- grep(paste0(prefix_master[j],".score"), colnames(chrom_file))
			substate_label <- state_score[[j]]$label[chrom_file[,pos_label]]
			substate_score <- state_score[[j]]$score[chrom_file[,pos_score]]
			chrom_comp[[j]]  <- cbind(substate_label,substate_score)
			colnames(chrom_comp[[j]]) <- colnames(chrom_file)[c(pos_label,pos_score)]
		}
		hmm_comp[[i]] <- do.call(cbind, chrom_comp)
	}
	hmm_comp <- do.call(rbind, hmm_comp)
	return(hmm_comp)
}

subcompartment <- function(data, diffdir, subnum) {

	diffdir <- paste0("DifferentialResult/",diffdir)
	prefix_master <- unique(data$prefix.master)
	if (file.exists(paste0(diffdir,"/fdr_result/differential.intra_sample_group.pcQnm.bedGraph"))) {
		compartment_file <- read.table(paste0(diffdir,"/fdr_result/differential.intra_sample_group.pcQnm.bedGraph"), h=T, as.is=T)
		hmm_comp <- hmmsegment(compartment_file, prefix_master, subnum)
		hmm_comp <- data.frame(compartment_file[,c(1:3,(ncol(compartment_file)-1))],hmm_comp)
		write.table(hmm_comp, file=paste0(diffdir,"/fdr_result/intra_sample_group.subcompartments.bedGraph"), row.names=F, col.names=T, sep="\t", quote=F)

		hmm_comp <- hmm_comp[hmm_comp$padj < 0.05, grep("state",colnames(hmm_comp))]
		mat_list <- list()
		for(i in 1:(length(prefix_master)-1)) {
			mat_list[[i]] <- table(hmm_comp[,i],hmm_comp[,(i+1)])
			rownames(mat_list[[i]]) <- paste0(prefix_master[i],".",rownames(mat_list[[i]]))
			colnames(mat_list[[i]]) <- paste0(prefix_master[(i+1)],".",colnames(mat_list[[i]]))
			mat_list[[i]] <- data.frame(mat_list[[i]])
		}
		mat_list <- data.frame(do.call(rbind, mat_list), stringsAsFactors = F)
		colnames(mat_list) <- c("source","target","value")
		nodes <- data.frame(name=unique(c(as.character(mat_list$source), as.character(mat_list$target))), stringsAsFactors = F)
		nodes[,"id"] <- 0:(nrow(nodes)-1)
		nodes[,"subcomp"] <- do.call(rbind,strsplit(as.character(nodes$name),"[.]"))[,-1]

		source <- list()
		target <- list()
		for(i in 1:nrow(mat_list)) {
			source[[i]] <- nodes[as.character(nodes$name) == as.character(mat_list$source[i]),]$id
			target[[i]] <- nodes[as.character(nodes$name) == as.character(mat_list$target[i]),]$id
		}

		mat_list[,"source"] <- unlist(source)
		mat_list[,"target"] <- unlist(target)
		mat_list <- mat_list[mat_list$value > quantile(mat_list$value,0),]
		network_df <- list(nodes=nodes,links=mat_list)
		ColourScal ='d3.scaleOrdinal() .range(["#D192D6","#A224AD","#F31D64","#FE433C","#B59CD9","#6A38B3","#3C50B1","#0095EF"])'
		sn <- networkD3::sankeyNetwork(Links = network_df$links, Nodes = network_df$nodes, Source = "source", Target = "target", Value = "value", NodeID = "subcomp", fontSize = 12, nodeWidth = 30, colourScale=ColourScal, sinksRight=T)
		jsobj <- htmlwidgets::JS('function(el, x) {','var cols_x = this.sankey.nodes().map(d => d.x).filter((v, i, a) => a.indexOf(v) === i);','var labels = [',paste0("\"",prefix_master,"\"",collapse=","),']', 'cols_x.forEach((d, i) => {', 'd3.select(el).select("svg")', '.append("text")', '.attr("x", d)', '.attr("y", 12)', '.text(labels[i]);', '})}')
		setwd(paste0(diffdir,"/fdr_result"))
		sn_header <- htmlwidgets::onRender (sn, jsCode=jsobj)
		networkD3::saveNetwork(sn_header, file = paste0("intra_sample_group.diff.subcompartments.transitions.html"))
		setwd("../../../")
	}

	if (file.exists(paste0(diffdir,"/fdr_result/differential.inter_sample_group.pcQnm.bedGraph"))) {
		compartment_file <- read.table(paste0(diffdir,"/fdr_result/differential.inter_sample_group.pcQnm.bedGraph"), h=T, as.is=T)
		hmm_comp <- hmmsegment(compartment_file, prefix_master, subnum)
		hmm_comp <- data.frame(compartment_file[,c(1:3,(ncol(compartment_file)-1))],hmm_comp, stringsAsFactors = F)
		write.table(hmm_comp, file=paste0(diffdir,"/fdr_result/inter_sample_group.subcompartments.bedGraph"), row.names=F, col.names=T, sep="\t", quote=F)
	
		hmm_comp <- hmm_comp[hmm_comp$padj < 0.05,grep("state",colnames(hmm_comp))]
		print (head(hmm_comp))
		mat_list <- list()
		for(i in 1:(length(prefix_master)-1)) {
			mat_list[[i]] <- table(hmm_comp[,i],hmm_comp[,(i+1)])
			rownames(mat_list[[i]]) <- paste0(prefix_master[i],".",rownames(mat_list[[i]]))
			colnames(mat_list[[i]]) <- paste0(prefix_master[(i+1)],".",colnames(mat_list[[i]]))
			mat_list[[i]] <- data.frame(mat_list[[i]])
		}
		mat_list <- data.frame(do.call(rbind, mat_list), stringsAsFactors = F)
		colnames(mat_list) <- c("source","target","value")
		print (head(mat_list))
		nodes <- data.frame(name=unique(c(as.character(mat_list$source), as.character(mat_list$target))), stringsAsFactors = F)
		nodes[,"id"] <- 0:(nrow(nodes)-1)
		nodes[,"subcomp"] <- do.call(rbind,strsplit(as.character(nodes$name),"[.]"))[,-1]

		source <- list()
		target <- list()
		for(i in 1:nrow(mat_list)) {
			source[[i]] <- nodes[as.character(nodes$name) == as.character(mat_list$source[i]),]$id
			target[[i]] <- nodes[as.character(nodes$name) == as.character(mat_list$target[i]),]$id
		}

		mat_list[,"source"] <- unlist(source)
		mat_list[,"target"] <- unlist(target)
		mat_list <- mat_list[mat_list$value > quantile(mat_list$value,0),]
		network_df <- list(nodes=nodes,links=mat_list)
		print (head(network_df))

		ColourScal ='d3.scaleOrdinal() .range(["#D192D6","#A224AD","#F31D64","#FE433C","#B59CD9","#6A38B3","#3C50B1","#0095EF"])'
		sn <- networkD3::sankeyNetwork(Links = network_df$links, Nodes = network_df$nodes, Source = "source", Target = "target", Value = "value", NodeID = "subcomp", fontSize = 12, nodeWidth = 30, colourScale=ColourScal, sinksRight=T)
		jsobj <- htmlwidgets::JS('function(el, x) {','var cols_x = this.sankey.nodes().map(d => d.x).filter((v, i, a) => a.indexOf(v) === i);','var labels = [',paste0("\"",prefix_master,"\"",collapse=","),']', 'cols_x.forEach((d, i) => {', 'd3.select(el).select("svg")', '.append("text")', '.attr("x", d)', '.attr("y", 12)', '.text(labels[i]);', '})}')
		setwd(paste0(diffdir,"/fdr_result"))
		sn_header <- htmlwidgets::onRender(sn, jsCode=jsobj)
		networkD3::saveNetwork(sn_header, file = paste0("inter_sample_group.diff.subcompartments.transitions.html"))
		setwd("../../../")
	}
}

### END ###


########################### Fithic calling ###########################

callfithiC <- function(x, fithic_path, python_path, resolution, df, diff_dir, u, dir_ovwt) {

	prefix <- as.character(df$prefix[x])
	folder <- paste0(diff_dir,"/fithic_run/",prefix,"_fithic")
	cmd <- paste0(python_path," ",fithic_path," -i ",folder,"/interactions.txt.gz -f ",folder,"/fragments.txt.gz -U ",as.integer(u)," -o ",folder,"/fithic_result -r ",as.integer(resolution))
	if (!file.exists(paste0(folder,"/fithic_result/FitHiC.spline_pass1.res",as.integer(resolution),".significances.txt.gz")) | dir_ovwt == TRUE) {
		cat (cmd,"\n")
		system(cmd, wait=T)
	} else {
		cat(paste0(folder,"/fithic_result/FitHiC.spline_pass1.res",as.integer(resolution),".significances.txt.gz file exists!\n"))	
	}
}


Rcpp::cppFunction('NumericVector MGLCC(std::vector<double> x, std::vector<double> y, int s) {
  double vec_array[s] = {};
  for(int i = 0; i < x.size(); i++) {
    int p = x[i] - 1;
    vec_array[p] = vec_array[p] + y[i];
  }
  NumericVector V = NumericVector(vec_array, vec_array+sizeof(vec_array)/sizeof(*vec_array));
  return V;
}')

formatconversion <- function(ijk_path, bed_path, prefix, diffdir, downsample=NA) {

	folder <- paste0(diffdir,"/fithic_run/",prefix,"_fithic")
	ifelse(!dir.exists(folder), dir.create(folder), print("folder exists"))
	if (!file.exists(paste0(folder,"/fragments.txt.gz"))) {
		cat ("Creating inputs for fithic run ",prefix,"\n")
		cat ("Reading the bed and interaction matrix file\n")
		bed_df <- read.table(bed_path, h=F)
		bed_df <- bed_df[,1:4]
		ijk_df <- as.data.frame(suppressWarnings(data.table::fread(ijk_path, h=F)))
		colnames(bed_df) <- c("chr","start","end","index")
		bed_df[,"extraField"] <- 0
		bed_df[,"mappable"] <- 1
	
		colnames(ijk_df) <- c("A","B","C")
		ijk_df[,"chr1"] <- NA
		ijk_df[,"fragmentMid1"] <- NA
		ijk_df[,"chr2"] <- NA
		ijk_df[,"fragmentMid2"] <- NA
	
		ratio <- 1
		if ((sum(ijk_df$C)/downsample) > 1 & !is.na(downsample)) {
			ratio = downsample/sum(ijk_df$C)
		}
		ijk_df[,"contactCount"] <- ceiling(c(ijk_df$C * ratio))
		ijk_df[ijk_df$contactCount > 0,]
	
		bed_df[,"mid"] <- bed_df$start + (bed_df$end[1]/2)
		bed_df$mid <- as.integer(bed_df$mid)
		mid_hash <- hashmap::hashmap(bed_df$index, bed_df$mid)
		chr_hash <- hashmap::hashmap(bed_df$index, as.character(bed_df$chr))
	
		ijk_df[,"chr1"] <- chr_hash[[ijk_df$A]]
		ijk_df[,"fragmentMid1"] <- mid_hash[[ijk_df$A]]
		ijk_df[,"chr2"] <- chr_hash[[ijk_df$B]]
		ijk_df[,"fragmentMid2"] <- mid_hash[[ijk_df$B]]
		ijk_df <- ijk_df[(ijk_df$chr1 == ijk_df$chr2),]
		ijk_df <- ijk_df[abs(ijk_df$fragmentMid2 - ijk_df$fragmentMid1) <= 10e6,]
		ijk_df <- na.omit(ijk_df)
		
		bed_df[,"correct_index"] <- 1:nrow(bed_df)
		index_hash <- hashmap::hashmap(bed_df$index, bed_df$correct_index)
		ijk_df[,"correct_A"] <- index_hash[[ijk_df$A]]
		ijk_df[,"correct_B"] <- index_hash[[ijk_df$B]]
		ijk_df$A <- ijk_df$A - (ijk_df$A - ijk_df$correct_A)
		ijk_df$B <- ijk_df$B - (ijk_df$B - ijk_df$correct_B)
		print (head(ijk_df))
		data.table::fwrite(ijk_df[,c(4:8)], file=paste0(folder,"/interactions.txt"), col.names=F, row.names=F, sep="\t", quote=F)
		R.utils::gzip(paste0(folder,"/interactions.txt"), destname=paste0(folder,"/interactions.txt.gz"), overwrite=T)

		self_cnt<- ijk_df[ijk_df$A == ijk_df$B,]
		ijk_df  <- ijk_df[ijk_df$A != ijk_df$B,]
		#self_cnt<- na.omit(self_cnt)
		#ijk_df  <- na.omit(ijk_df) 
		pos_vec <- c(ijk_df$A,ijk_df$B,self_cnt$A)
		cnt_vec <- c(ijk_df$C,ijk_df$C,self_cnt$C)
		
		cat ("Started calculating Marginalized Contact Count\n")
		index_arr <- MGLCC(pos_vec, cnt_vec, ifelse(nrow(bed_df) < max(pos_vec), max(pos_vec), nrow(bed_df)))
		bed_df[,"marginalizedContactCount"] <- index_arr
		print (head(bed_df))
		data.table::fwrite(bed_df[,c("chr","extraField","mid","marginalizedContactCount","mappable")], file=paste0(folder,"/fragments.txt"), col.names=F, row.names=F, sep="\t", quote=F)
		R.utils::gzip(paste0(folder,"/fragments.txt"), destname=paste0(folder,"/fragments.txt.gz"), overwrite=T)
		rm(ijk_df)
		rm(mid_hash)
		rm(chr_hash)
	} else {
		cat ("Fithic file already exists for ",prefix,", skipping\n")		
	}
}


fithicformat <- function(data, diffdir, fithicpath, pythonpath, fdr_thr, dist_thr, dirovwt, cl_sthread) {

	diffdir <- paste0("DifferentialResult/",diffdir)
	if (!dir.exists(paste0(diffdir,"/fithic_run"))) {
		dir.create(paste0(diffdir,"/fithic_run"))
	}

	for(i in 1:nrow(data)) {
		formatconversion(normalizePath(data$mat[i]), normalizePath(data$bed[i]), as.character(data$prefix[i]), diffdir)
	}
	
	resolution <- read.table(data$bed[1], h=F, as.is=T)
	resolution <- resolution$V3[1] - resolution$V2[1]
	resolution <- as.integer(resolution)

	if (!is.na(cl_sthread)) {
		parallel::clusterExport(cl_sthread, "callfithiC")
		invisible(parallel::parLapply(cl_sthread, c(1:nrow(data)), callfithiC, fithicpath, pythonpath, resolution, data, diffdir, dist_thr, dirovwt))
	} else {
		invisible(lapply(c(1:nrow(data)), callfithiC, fithicpath, pythonpath, resolution, data, diffdir, dist_thr, dirovwt))
	}

	prefix_master <- unique(as.character(data$prefix.master))
	mat_sample <- list()
	ids_sample <- list()
	for(i in 1:length(prefix_master)) {
		print (i)
		data_rep <- data[data$prefix.master==as.character(prefix_master[i]),]
		mat_rep  <- list()
		ids_rep  <- list()
		for(j in 1:nrow(data_rep)) {
			mat_rep[[j]] <- data.table::fread(paste0("gzip -dc ",diffdir,"/fithic_run/",data_rep$prefix[j],"_fithic/fithic_result/FitHiC.spline_pass1.res",as.integer(resolution),".significances.txt.gz"), h=T)
			colnames(mat_rep[[j]])[c(6,7)] <- c("pval","qval")
			mat_rep[[j]] <- mat_rep[[j]][mat_rep[[j]]$qval < fdr_thr,]
			mat_rep[[j]][,"id"]  <- paste0(mat_rep[[j]]$chr1,"_",as.integer(mat_rep[[j]]$fragmentMid1-(resolution/2)),"_",mat_rep[[j]]$chr2,"_",as.integer(mat_rep[[j]]$fragmentMid2-(resolution/2)))
			mat_rep[[j]][,"sig"] <- 1
			data.table::setkey(mat_rep[[j]],id)
			ids_rep[[j]] <- as.character(unlist(mat_rep[[j]][,11]))
		}
		ids_rep <- sort(unique(as.character(unlist(ids_rep))))
		mat_sample[[i]] <- matrix(0, length(ids_rep), nrow(data_rep))
		rownames(mat_sample[[i]]) <- ids_rep
		colnames(mat_sample[[i]]) <- as.character(data_rep$prefix)
		for(j in 1:nrow(data_rep)) {
			mat_sample[[i]][,j] <- as.integer(unlist(mat_rep[[j]][.(ids_rep)][,12]))
		}
		rm(mat_rep)
		mat_sample[[i]][is.na(mat_sample[[i]])] <- 0
		freq_thr <- ifelse(nrow(data_rep) == 1, 0, 0.5)
		mat_sample[[i]] <- mat_sample[[i]][which((apply(mat_sample[[i]], 1, sum)/nrow(data_rep)) >= freq_thr),]
		if (nrow(data_rep) > 1) {
			ids_sample[[i]] <- rownames(mat_sample[[i]])
		} else {
			ids_sample[[i]] <- names(mat_sample[[i]])
		}
		mat_sample[[i]] <- data.table::as.data.table(mat_sample[[i]])
		mat_sample[[i]][,"id"]  <- ids_sample[[i]]
		mat_sample[[i]][,"sig"] <- 1
		data.table::setkey(mat_sample[[i]],id)
	}
	ids_sample <- sort(unique(as.character(unlist(ids_sample))))
	mat <- matrix(0, length(ids_sample), length(prefix_master))
	rownames(mat) <- ids_sample
	colnames(mat) <- prefix_master
	for(i in 1:length(prefix_master)) {
		mat[,i] <- as.integer(unlist(mat_sample[[i]][.(ids_sample)]$sig))
	}
	mat[is.na(mat)] <- 0
	write.table(mat, file=paste0(diffdir,"/fithic_run/FithicResult.txt"), row.names=T, col.names=T, sep="\t", quote=F)
}

### END ###

########################### Differential Loops ###########################

compartmentLoop <- function(data, diffdir, fdr.thr, count.thr=5, loop.dist=2e6) {

	prefix_master <- unique(data$prefix.master)
	diff_res   <- read.table(paste0("DifferentialResult/",diffdir,"/fdr_result/differential.intra_sample_group.Filtered.pcQnm.bedGraph"), h=T, as.is=T)
	diff_hash  <- hashmap::hashmap(as.character(paste0(diff_res$chr,"_",diff_res$start)), rep(1, nrow(diff_res)))
	chr_list   <- unique(diff_res$chr)
	resolution <- read.table(data$bed[1], h=F, as.is=T)
	resolution <- resolution$V3[1] - resolution$V2[1]
	
	mat_diff_loop <- list()
	for(i in 1:length(chr_list)) {
		region_list <- list()
		diff_sample_OE <- list()
		diff_sample_Count <- list()
		for(j in 1:length(prefix_master)) {
			data_rep <- data[data$prefix.master==prefix_master[j],]
			diff_rep <- list()
			n <- 1
			for(k in 1:nrow(data_rep)) {
				cat ("Reading ",paste0(data_rep$prefix[k],"_pca/intra_pca/",data_rep$prefix[k],"_mat/",chr_list[i],".txt"),"\n")
				chr_data <- data.table::fread(paste0(data_rep$prefix[k],"_pca/intra_pca/",data_rep$prefix[k],"_mat/",chr_list[i],".txt"), h=T)
				chr_data <- chr_data[abs(chr_data$pos2 - chr_data$pos1) <= loop.dist,]
				chr_data[,"match1"] <- diff_hash[[paste0(chr_data$chr1,"_",chr_data$pos1)]]
				chr_data[,"match2"] <- diff_hash[[paste0(chr_data$chr2,"_",chr_data$pos2)]]
				chr_data$match1[is.na(chr_data$match1)] <- 0
				chr_data$match2[is.na(chr_data$match2)] <- 0
				diff_rep[[n]] <- chr_data[(chr_data$match1==1 | chr_data$match2==1),]
				n <- n + 1
			}
			diff_rep <- data.table::rbindlist(diff_rep)
			diff_sample_OE[[j]] <- aggregate(WeightOE ~ chr1 + pos1 + chr2 + pos2, mean, data=diff_rep)
			diff_sample_OE[[j]][,"id"] <- paste0(diff_sample_OE[[j]]$chr1,"_",diff_sample_OE[[j]]$pos1,"_",diff_sample_OE[[j]]$chr2,"_",diff_sample_OE[[j]]$pos2)
			diff_sample_Count[[j]] <- aggregate(Weight ~ chr1 + pos1 + chr2 + pos2, mean, data=diff_rep)
			diff_sample_Count[[j]][,"id"] <- paste0(diff_sample_Count[[j]]$chr1,"_",diff_sample_Count[[j]]$pos1,"_",diff_sample_Count[[j]]$chr2,"_",diff_sample_Count[[j]]$pos2)
			region_list[[j]] <- paste0(diff_sample_OE[[j]]$chr1,"_",diff_sample_OE[[j]]$pos1,"_",diff_sample_OE[[j]]$chr2,"_",diff_sample_OE[[j]]$pos2)
			diff_sample_OE[[j]] <- data.table::as.data.table(diff_sample_OE[[j]])
			diff_sample_Count[[j]] <- data.table::as.data.table(diff_sample_Count[[j]])
			data.table::setkey(diff_sample_OE[[j]], "id")
			data.table::setkey(diff_sample_Count[[j]], "id")
		}

		region_list   <- unique(unlist(region_list))
		mat <- matrix(0, length(region_list), (2 * length(prefix_master)))
		colnames(mat) <- c(paste0(prefix_master,".WeightOE"), paste0(prefix_master,".Count"))
		rownames(mat) <- region_list
		for(n in 1:length(diff_sample_OE)) {
			mat[,n] <- as.numeric(unlist(diff_sample_OE[[n]][.(as.character(region_list))][,5]))
			mat[,(length(prefix_master) + n)] <- as.numeric(unlist(diff_sample_Count[[n]][.(as.character(region_list))][,5]))
		}

		mat[is.na(mat)] <- 0
		mat <- as.data.frame(mat)
		mat <- mat[which(apply(mat[,grep("Count",colnames(mat))], 1, min) > count.thr),]
		mat <- mat[which(apply(mat[,grep("WeightOE",colnames(mat))], 1, min) < 10),]
		mat <- mat[which(apply(mat[,grep("WeightOE",colnames(mat))], 1, max) > 0.5),]
		mat[,"OE.mean"] <- apply(mat[,grep("WeightOE", colnames(mat))], 1, mean)
    	mat.cen <- round(log2(mat[,grep("WeightOE", colnames(mat))]/mat$OE.mean), 5)
    	mat[,"maha"] <- round(mahalanobis(mat.cen, c(rep(0, ncol(mat.cen))), robust::covRob(mat.cen)$cov),5)
    	mat_diff_loop[[i]] <- mat
    	mat_diff_loop[[i]][,"id"] <- rownames(mat)
	}

	mat_diff_loop <- as.data.frame(data.table::rbindlist(mat_diff_loop))
	mat_diff_loop[,"padj"] <- p.adjust(pchisq(mat_diff_loop$maha, df=length(prefix_master)-1, lower.tail=F),"BH")
	mat_diff_loop[,c("chr1","start1","chr2","start2")] <- as.data.frame(do.call(rbind,strsplit(as.character(mat_diff_loop$id),"_")), stringsAsFactors=F)
	mat_diff_loop <- mat_diff_loop[mat_diff_loop$start1 != mat_diff_loop$start2,]
	write.table(mat_diff_loop, file=paste0("DifferentialResult/",diffdir,"/fdr_result/differential.intra_compartmentLoops.txt"), row.names=F, col.names=T, sep="\t", quote=F)

	if (file.exists(paste0("DifferentialResult/",diffdir,"/fithic_run/FithicResult.txt"))) {
		
		cat ("Found fithic significant loop calling result. Using it to futher filter differential loops\n")
		fithicresult  <- read.table(paste0("DifferentialResult/",diffdir,"/fithic_run/FithicResult.txt"), h=T, as.is=T)
		fithicresult  <- fithicresult[apply(fithicresult, 1, sum) < ncol(fithicresult),]
		mat_diff_loop <- data.table::as.data.table(mat_diff_loop)
		data.table::setkey(mat_diff_loop, id)
		mat_diff_loop <- na.omit(mat_diff_loop[.(rownames(fithicresult))])
		write.table(cbind(mat_diff_loop,fithicresult[mat_diff_loop$id,]), file=paste0("DifferentialResult/",diffdir,"/fdr_result/differential.intra_compartmentLoops.FithicSignificant.txt"), row.names=F, col.names=T, sep="\t", quote=F)

		mat_diff_loop <- mat_diff_loop[mat_diff_loop$padj < fdr.thr,]
		mat_diff_loop <- as.data.frame(mat_diff_loop)
		mat_diff_loop[,"end1"] <- as.integer(mat_diff_loop$start1) + as.integer(resolution)
		mat_diff_loop[,"end2"] <- as.integer(mat_diff_loop$start2) + as.integer(resolution)
		mat_diff_loop[,"logPadj"] <- -log10(mat_diff_loop$padj)
		mat_diff_loop$logPadj[mat_diff_loop$logPadj == "Inf"] <- max(mat_diff_loop$logPadj[mat_diff_loop$logPadj != "Inf"]) + 1
		print (head(mat_diff_loop))
		mat_diff_loop <- mat_diff_loop[,c("chr1","start1","end1","chr2","start2","end2","logPadj")]
		write.table(mat_diff_loop, file=paste0("DifferentialResult/",diffdir,"/fdr_result/differential.intra_compartmentLoops.bedpe"), row.names=F, col.names=T, sep="\t", quote=F)

	} else {

		cat ("FithicResult.txt file is missing. Using distance normalized Observed by Expected ratio to further filter the differential loops\n")
		mat_diff_loop <- mat_diff_loop[mat_diff_loop$padj < fdr.thr,]
		mat_diff_loop <- mat_diff_loop[which(apply(mat_diff_loop[,grep("WeightOE",colnames(mat_diff_loop))], 1, max) > 2),]
		mat_diff_loop[,"end1"] <- as.integer(mat_diff_loop$start1) + as.integer(resolution)
		mat_diff_loop[,"end2"] <- as.integer(mat_diff_loop$start2) + as.integer(resolution)
		mat_diff_loop[,"logPadj"] <- -log10(mat_diff_loop$padj)
		mat_diff_loop$logPadj[mat_diff_loop$logPadj == "Inf"] <- max(mat_diff_loop$logPadj[mat_diff_loop$logPadj != "Inf"]) + 1
		print (head(mat_diff_loop))
		mat_diff_loop <- mat_diff_loop[,c("chr1","start1","end1","chr2","start2","end2","logPadj")]
		write.table(mat_diff_loop, file=paste0("DifferentialResult/",diffdir,"/fdr_result/differential.intra_compartmentLoops.bedpe"), row.names=F, col.names=T, sep="\t", quote=F)

	}
}

### END ###

########################### Visualization ###########################

htmlheader <- function() {
cat ("<!DOCTYPE html>
<html lang=\"en\">
<head>
    <meta charset=\"utf-8\">
    <meta content=\"IE=edge\" http-equiv=\"X-UA-Compatible\">
    <meta content=\"width=device-width, initial-scale=1, maximum-scale=1, user-scalable=no\" name=\"viewport\">
    <meta content=\"\" name=\"description\">
    <meta content=\"\" name=\"",paste0(Sys.info()["user"][[1]]),"\">
    <link href=https://igv.org/web/release/2.9.1/examples/img/favicon.ico rel=\"igv icon\">
    <title>IGV - Dev</title>
    <script src=\"https://cdn.jsdelivr.net/npm/igv@2.9.1/dist/igv.min.js\"></script>
</head>

<body>
<h1><u>dcHiC: Differential Compartment calling from Hi-C data</u></h1>
<p style=\"color: #C86400\"><b>A compartment</b></p>
<p style=\"color: #0064C8\"><b>B compartment</b></p>
<p style=\"color: #99FFCC\"><b>Mdist: Mahalanobis distance score to represent outlierness of the bin</b></p>
<p style=\"color: #009900\"><b>log10Padj: P adjusted value of corresponding Mdist score</b></p>
<p style=\"color: #994C00\"><b>dZsc: Distance Zscore that represent biological variability as compared to technical variability</b></p>
<div id=\"igvDiv\" style=\"padding-top: 10px;padding-bottom: 10px; border:1px solid lightgray\"></div>\n")
}

scriptbody <- function(genome,start=F,end=F) {

  if (start) {
  cat ("<script type=\"text/javascript\">
    document.addEventListener(\"DOMContentLoaded\", function () {
      var igvDiv = document.getElementById(\"igvDiv\");
      var options = {
        locus: '19:49301000-49305700',
        genome:",paste0("\"",genome,"\""),",
	tracks: [\n")
  }
  if (end) {
  cat ("
       ]
     }; igv.createBrowser(igvDiv, options)
                .then(function (browser) {
                    console.log(\"Created IGV browser\")});
   })
  </script>
 </body>
</html>\n")
  }
}


htmlbody <- function(file, folder, genome, fdr_thr, numberclust) {

 	if (!dir.exists(paste0(folder,"/data"))) {
		dir.create(paste0(folder,"/data"))
	}
	scriptbody(genome,start=T)
  
	for(i in 1:nrow(file)) {
    	if (file$group[i] == "compartment") {

      		compbdg <- read.table(as.character(file$file[i]), h=T, as.is=T, comment.char = "")
      		colnames(compbdg)[1] <- "chr"

      		color = "00,100,0"
      		altcolor = "0,100,200"
      		if (ncol(file) == 4) {
	       		color = strsplit(as.character(file$color[i]),":")[[1]][1]
	       		altcolor = strsplit(as.character(file$color[i]),":")[[1]][2]
      		}

      		## Write the PC coordinates
      		for(j in 4:(ncol(compbdg)-3)) {
			cat("# locus chr19:49302001-49304701\n# refGene encodeRegions\n# zero-based, half-open coords
track name=\"",colnames(compbdg)[j]," PC\" description=\"BedGraph format\" visibility=full color=",color," altColor=",altcolor," priority=20 plotType=\"points\"\n", 
			file=paste0(folder,"/data/",colnames(compbdg)[j],".PC.bedGraph"))

      		write.table(compbdg[,c(1:3,j)], file=paste0(folder,"/data/",colnames(compbdg)[j],".PC.bedGraph"), row.names=F, col.names=F,sep="\t", quote=F, append=T)
      		R.utils::gzip(paste0(folder,"/data/",colnames(compbdg)[j],".PC.bedGraph"), paste0(folder,"/data/",colnames(compbdg)[j],".PC.bedGraph.gz"), overwrite=T)
      		cmd <- paste0("create_datauri ",folder,"/data/",colnames(compbdg)[j],".PC.bedGraph.gz")
      		datauri <- system(cmd, wait=T, intern=T)
        	cat ("
             	{
                	name: ",paste0("\'",colnames(compbdg)[j],".PC',"),"
                	url: ",paste0("\"",datauri,"\""),",
                	indexed: false,
                	format: \"bedGraph\"
             	},\n")       
      		}

 			trackmaster <- paste0(colnames(compbdg)[4:(ncol(compbdg)-1)],collapse="_vs_") 
     		## Write Padjusted values
			cat("# locus chr19:49302001-49304701\n# refGene encodeRegions\n# zero-based, half-open coords
track name=\"log10Padj\" description=\"BedGraph format\" visibility=full color=0,153,0 priority=20 plotType=\"points\"\n",
			file=paste0(folder,"/data/differential_compartment.log10Padj.bedGraph"))    

			compbdg_tmp <- compbdg[compbdg$padj > -log10(fdr_thr) & compbdg$dist_clust > numberclust,]
     		write.table(compbdg_tmp[,c(1:3,(ncol(compbdg_tmp)-1))], file=paste0(folder,"/data/differential_compartment.log10Padj.bedGraph"), row.names=F, col.names=F,sep="\t", quote=F, append=T)
     		R.utils::gzip(paste0(folder,"/data/differential_compartment.log10Padj.bedGraph"), paste0(folder,"/data/differential_compartment.log10Padj.bedGraph.gz"), overwrite=T)
     		cmd <- paste0("create_datauri ",folder,"/data/differential_compartment.log10Padj.bedGraph.gz")
     		datauri <- system(cmd, wait=T, intern=T)
        	cat ("
             	{
                	name: ",paste0("\'log10Padj\',"),"
                	url: ",paste0("\"",datauri,"\""),",
                	indexed: false,
                	format: \"bedGraph\"
             	},\n")

      		## Write Replicate weight
			#cat("# locus chr19:49302001-49304701\n# refGene encodeRegions\n# zero-based, half-open coords
#track name=\"Replicate weight\" description=\"BedGraph format\" visibility=full color=153,148,194 priority=20 plotType=\"points\"\n",
			#file=paste0(folder,"/data/",trackmaster,".ReplicateWeight.bedGraph"))    

     		#write.table(compbdg[,c(1:3,(ncol(compbdg)-3))], file=paste0(folder,"/data/",trackmaster,".ReplicateWeight.bedGraph"), row.names=F, col.names=F,sep="\t", quote=F, append=T)
     		#R.utils::gzip(paste0(folder,"/data/",trackmaster,".ReplicateWeight.bedGraph"), paste0(folder,"/data/",trackmaster,".ReplicateWeight.bedGraph.gz"), overwrite=T)
     		#cmd <- paste0("create_datauri ",folder,"/data/",trackmaster,".ReplicateWeight.bedGraph.gz")
     		#datauri <- system(cmd, wait=T, intern=T)
        	#cat ("
            # 	{
            #    	name: ",paste0("\'Replicate weight\',"),"
            #    	url: ",paste0("\"",datauri,"\""),",
            #    	indexed: false,
            #    	format: \"bedGraph\"
            # 	},\n")


      		## Write Sample mahalanobis
			cat("# locus chr19:49302001-49304701\n# refGene encodeRegions\n# zero-based, half-open coords
track name=\"Mahalanobis distance\" description=\"BedGraph format\" visibility=full color=284,102,84 priority=20 plotType=\"points\"\n",
			file=paste0(folder,"/data/differential_compartment.Mahalanobis.bedGraph"))    

     		write.table(compbdg[,c(1:3,(ncol(compbdg)-2))], file=paste0(folder,"/data/differential_compartment.Mahalanobis.bedGraph"), row.names=F, col.names=F,sep="\t", quote=F, append=T)
     		R.utils::gzip(paste0(folder,"/data/differential_compartment.Mahalanobis.bedGraph"), paste0(folder,"/data/differential_compartment.Mahalanobis.bedGraph.gz"), overwrite=T)
     		cmd <- paste0("create_datauri ",folder,"/data/differential_compartment.Mahalanobis.bedGraph.gz")
     		datauri <- system(cmd, wait=T, intern=T)
        	cat ("
             	{
                	name: ",paste0("\'Mahalanobis distance\',"),"
                	url: ",paste0("\"",datauri,"\""),",
                	indexed: false,
                	format: \"bedGraph\"
             	}")


        	## Write Glosh score
        	cat("# locus chr19:49302001-49304701\n# refGene encodeRegions\n# zero-based, half-open coords
track name=\"Glosh Score\" description=\"BedGraph format\" visibility=full color=184,162,84 priority=20 plotType=\"points\"\n",
			file=paste0(folder,"/data/differential_compartment.Mahalanobis.bedGraph"))    

     		#write.table(compbdg[,c(1:3,(ncol(compbdg)))], file=paste0(folder,"/data/differential_compartment.GloshScore.bedGraph"), row.names=F, col.names=F,sep="\t", quote=F, append=T)
     		#R.utils::gzip(paste0(folder,"/data/differential_compartment.GloshScore.bedGraph"), paste0(folder,"/data/differential_compartment.GloshScore.bedGraph.gz"), overwrite=T)
     		#cmd <- paste0("create_datauri ",folder,"/data/differential_compartment.GloshScore.bedGraph.gz")
     		#datauri <- system(cmd, wait=T, intern=T)
        	#cat ("
            # 	{
            #    	name: ",paste0("\'Glosh Score\',"),"
            #    	url: ",paste0("\"",datauri,"\""),",
            #    	indexed: false,
            #    	format: \"bedGraph\"
            # 	}")

      		if (i < nrow(file)) { 
       			cat (",\n") 
      		}
    	}	else if (file$group[i] == "seg") {

    		if (file$name[i] != "Segment_Additional") {
				R.utils::gzip(as.character(file$file[i]), paste0(folder,"/data/intra_",file$name[i],".seg.gz"), overwrite=T)
				cmd <- paste0("create_datauri ",folder,"/data/intra_",file$name[i],".seg.gz")
     			datauri <- system(cmd, wait=T, intern=T)
        		cat ("
            	 	{
            	    	name: ",paste0("\'",file$name[i],"\',"),"
            	    	url: ",paste0("\"",datauri,"\""),",
            	    	indexed: false,
            	    	isLog: true,
            	    	height: 25,
            	    	maxHeight: 50,
            	    	format: \"seg\"
            	 	},\n")
            } else if (file$name[i] == "Segment_Additional") {
            	#R.utils::gzip(as.character(file$file[i]), paste0(folder,"/data/",basename(as.character(file$file[i])),".seg.gz"), overwrite=T)
            	cmd <- paste0("cp ",as.character(file$file[i])," ",folder,"/data/",basename(as.character(file$file[i])))
            	system(cmd, wait=T, intern=T)
            	cmd <- paste0("create_datauri ",folder,"/data/",basename(as.character(file$file[i])))
            	datauri <- system(cmd, wait=T, intern=T)
        		cat ("
            	 	{
            	    	name: ",paste0("\'",basename(as.character(file$file[i])),"\',"),"
            	    	url: ",paste0("\"",datauri,"\""),",
            	    	indexed: false,
            	    	isLog: true,
            	    	height: 25,
            	    	maxHeight: 50,
            	    	format: \"seg\"
            	 	},\n")
            }

    	} 	else if (file$group[i] == "bedGraph") {

        	bdg <- read.table(as.character(file$file[i]),h=F)
        	trackname <- as.character(file$name[i])
			color = "255,102,178"
     		altcolor = "102,102,255"
      		if (ncol(file) == 4) {
          		color = strsplit(as.character(file$color[i]),":")[[1]][1]
          		if (strsplit(as.character(file$color[i]),":")[[1]][2] != "0,0,0") {
            		altcolor = strsplit(as.character(file$color[i]),":")[[1]][2]
          		}
      		}

		cat("# locus chr19:49302001-49304701\n# refGene encodeRegions\n# zero-based, half-open coords
track name=\"",trackname,"\" description=\"BedGraph format\" visibility=full color=",color," altColor=",altcolor," priority=20 plotType=\"points\"\n",
			file=paste0(folder,"/data/",trackname,".bedGraph"))
        	write.table(bdg, file=paste0(folder,"/data/",trackname,".bedGraph"),row.names=F, col.names=F,sep="\t", quote=F, append=T)
        	R.utils::gzip(paste0(folder,"/data/",trackname,".bedGraph"), paste0(folder,"/data/",trackname,".bedGraph.gz"), overwrite=T)
        	cmd <- paste0("create_datauri ",folder,"/data/",trackname,".bedGraph.gz")
        	datauri <- system(cmd, wait=T, intern=T)
        	cat ("
        	     {
            	    name: ",paste0("\'",trackname,"',"),"
            	    url: ",paste0("\"",datauri,"\""),",
            	    indexed: false,
            	    format: \"bedGraph\"
             	}")
        	if (i < nrow(file)) {
        	  cat (",\n")
        	}
    	} else if (file$group[i] == "bedInt") {

    		if (file$name[i] == "DifferentialLoops") {
    			cmd <- paste0("cp ",as.character(file$file[i])," ",folder,"/data/differential.intra_compartmentLoops.bedpe")
    			system(cmd, wait=T)
    			cmd <- paste0("create_datauri ",folder,"/data/differential.intra_compartmentLoops.bedpe")
    			datauri <- system(cmd, wait=T, intern=T)
    			cat ("
        		     {
            		    name: ",paste0("\'Differential Loops\',"),"
            		    url: ",paste0("\"",datauri,"\""),",
            		    type: \"interaction\",
            		    format: \"bedpe\"
            	 	}")

    			if (i < nrow(file)) {
        	 		cat (",\n")
        		}
    		} else if (file$name[i] == "DifferentialLoops_Additional") {
    			cmd <- paste0("cp ",as.character(file$file[i])," ",folder,"/data/",basename(as.character(file$file[i])))
    			system(cmd, wait=T)
    			cmd <- paste0("create_datauri ",folder,"/data/",basename(as.character(file$file[i])))
    			datauri <- system(cmd, wait=T, intern=T)
    			cat ("
        		     {
            		    name: ",paste0("\'",basename(as.character(file$file[i]))," Differential Loops\',"),"
            		    url: ",paste0("\"",datauri,"\""),",
            		    type: \"interaction\",
            		    format: \"bedpe\"
            	 	}")

    			if (i < nrow(file)) {
        	 		cat (",\n")
        		}
    		}
    	}
  	} 
  	scriptbody(genome,end=T) 
}

generateTrackfiles <- function(data, diffdir, genome, bdgfile, pcgrp="pcQnm", fdr_thr, numberclust) {

	diffdir <- paste0("DifferentialResult/",diffdir)
	if (!dir.exists(paste0(diffdir,"/viz"))) {
		dir.create(paste0(diffdir,"/viz"))	
	}
	if (!dir.exists(paste0(diffdir,"/viz/files"))) {
		dir.create(paste0(diffdir,"/viz/files"))
	}

	# Intra IGV html file generation
	prefix_master <- unique(data$prefix.master)
	if (file.exists(paste0(diffdir,"/fdr_result/differential.intra_sample_group.",pcgrp,".bedGraph"))) {
		compartment_file <- read.table(paste0(diffdir,"/fdr_result/differential.intra_sample_group.",pcgrp,".bedGraph"), h=T, as.is=T)
		compartment_file <- compartment_file[order(compartment_file$chr, compartment_file$start),c("chr","start","end",prefix_master,"sample_maha","padj","dist_clust")]
		print (head(compartment_file))
		compartment_file$padj[compartment_file$padj == 0] <- min(compartment_file$padj[compartment_file$padj > 0])
		compartment_file$padj <- -log10(compartment_file$padj)
		write.table(compartment_file[,1:(ncol(compartment_file))], file=paste0(diffdir,"/viz/files/intra_compartment.bedGraph"), row.names=F, col.names=T, sep="\t", quote=F)

		intra_pcori <- list.files(paste0(diffdir,"/",pcgrp), pattern="intra_sample", full.names=T)
		pcori_files <- list()
		for(i in 1:length(intra_pcori)) {
			pcori_files[[i]] <- read.table(intra_pcori[i], h=T, as.is=T)
		}

		pcori_files <- do.call(rbind, pcori_files)
		pcori_files <- pcori_files[order(pcori_files$chr, pcori_files$start),]
		for(i in 1:(ncol(pcori_files) - 3)) {
			write.table(pcori_files[,c(1:3,(i+3))], file=paste0(diffdir,"/viz/files/intra_",colnames(pcori_files)[i+3],"_PC.bedGraph"), row.names=F, col.names=F, sep="\t", quote=F)
		}

		if (file.exists(paste0(diffdir,"/fdr_result/intra_sample_group.subcompartments.bedGraph"))) {
			subcomp <- read.table(paste0(diffdir,"/fdr_result/intra_sample_group.subcompartments.bedGraph"), h=T, as.is=T)
			for(j in 1:length(prefix_master)) {
				subcomp_temp <- subcomp[,c("chr","start","end",paste0(prefix_master[j],".score"))]
				colnames(subcomp_temp)  <- c("Chromosome","Start","End","Segment_Mean")
				subcomp_temp[,"Sample"] <- prefix_master[j]
				subcomp_temp[,"Num_Probes"] <- 1
				subcomp_temp <- subcomp_temp[,c("Sample","Chromosome","Start","End","Num_Probes","Segment_Mean")]
				write.table(subcomp_temp, file=paste0(diffdir,"/viz/files/intra_",prefix_master[j],"_subcomp.seg"), row.names=F, col.names=T, sep="\t", quote=F)
			}
		}

		
		if (file.exists(paste0(diffdir,"/fdr_result/differential.intra_compartmentLoops.bedpe"))) {
			bedpe <- read.table(paste0(diffdir,"/fdr_result/differential.intra_compartmentLoops.bedpe"), h=T, as.is=T)
			bigInt<- data.frame(
					chrom=bedpe$chr1,
					chromStart=bedpe$start1,
					chromEnd=bedpe$end1,
					targetChrom=bedpe$chr2,
                    targetStart=bedpe$start2,
                    targetEnd=bedpe$end2,
					value=round(bedpe$logPadj,3)
				)

			bigInt <- bigInt[order(bigInt$chrom, bigInt$chromStart),]
			write.table(bigInt, file=paste0(diffdir,"/viz/files/differential.intra_compartmentLoops.bedpe"), row.names=F, col.names=F, sep="\t", quote=F)
		}	

		file <- list()
		file[[1]] <- data.frame(
					file=paste0("./files/intra_compartment.bedGraph"),
					name="Selected",
					group="compartment",
					color="234,100,0:122,16,180"
				)

		comp_colA <- t(col2rgb(colorspace::diverge_hcl(length(prefix_master), c=100, l=c(50,90), power=0.5)))
		comp_colB <- t(col2rgb(colorspace::heat_hcl(length(prefix_master), h=c(80,-100), l=c(75,40), c=c(40,80), power=0.5)))
		#k <- length(prefix_master)
		#l <- 2
		#for(i in 1:length(prefix_master)) {
		#	data_rep <- data[data$prefix.master== as.character(prefix_master[i]),] 
		#	print (data_rep)
		#	if (comp_colA[k,1] > comp_colB[k,1]) {
		#		red_col <- paste0(comp_colA[k,], collapse=",")
		#		blu_col <- paste0(comp_colB[k,], collapse=",")
		#	} else {
		#		red_col <- paste0(comp_colB[k,], collapse=",")
		#		blu_col <- paste0(comp_colA[k,], collapse=",")
		#	}
		#	for(j in 1:nrow(data_rep)) {
		#		file[[l]] <- data.frame(
		#					file=paste0("./files/intra_",data_rep$prefix[j],"_PC.bedGraph"),
		#					name=data_rep$prefix[j],
		#					group="bedGraph",
		#					color=paste0(red_col,":",blu_col)
		#				)
		#		l <- l + 1
		#	}
		#	k <- k - 1
		#}
		
		if (file.exists(paste0(diffdir,"/fdr_result/intra_sample_group.subcompartments.bedGraph"))) {
			for(j in 1:length(prefix_master)) {
				file[[length(file)+j]] <- data.frame(
							file=paste0("./files/intra_",prefix_master[j],"_subcomp.seg"),
							name=paste0(prefix_master[j],".subcomp"),
							#name="Segment",
							group="seg",
							color="0:0"
						)
			}
		}

		if (file.exists(paste0(diffdir,"/fdr_result/differential.intra_compartmentLoops.bedpe"))) {
			file[[length(file)+1]] <- data.frame(
						file=paste0("./files/differential.intra_compartmentLoops.bedpe"),
						name="DifferentialLoops",
						group="bedInt",
						color="0:0"
				)
		}


		if (!is.na(bdgfile)) {
			bdgfile <- read.table(bdgfile, h=F, as.is=T)
			colnames(bdgfile) <- c("file", "name", "type", "color")
			for(b in 1:nrow(bdgfile)) {
				if (bdgfile$type[b] == "bedGraph") {
					cmd <- paste0("cp ",bdgfile$file[b]," ",diffdir,"/viz/files/")
					system(cmd, wait=T)
					file[[length(file)+1]] <- data.frame(
							file=paste0("./files/",basename(bdgfile$file[b])),
							name=bdgfile$name[b],
							group="bedGraph",
							color=bdgfile$color[b]
					)
				} else if (bdgfile$type[b] == "bedInt") {
					cmd <- paste0("cp ",bdgfile$file[b]," ",diffdir,"/viz/files/")
					system(cmd, wait=T)
					file[[length(file)+1]] <- data.frame(
							file=paste0("./files/",basename(bdgfile$file[b])),
							name="DifferentialLoops_Additional",
							group="bedInt",
							color="0:0"
					)
				}
				else if (bdgfile$type[b] == "seg") {
					cmd <- paste0("cp ",bdgfile$file[b]," ",diffdir,"/viz/files/")
					system(cmd, wait=T)
					file[[length(file)+1]] <- data.frame(
							file=paste0("./files/",basename(bdgfile$file[b])),
							name="Segment_Additional",
							group="seg",
							color="0:0"
					)
				}
			}
		}

		file <- as.data.frame(do.call(rbind, file), stringsAsFactors=F)
		print (file)

		current_path <- getwd()
		setwd(paste0(diffdir,"/viz/"))
		folder <- "vizIGV_intra"
		if (!file.exists(folder)) {
 			dir.create(folder)
		}
		sink(paste0(folder,"/intra_igv.html"))
		htmlheader()
		htmlbody(file, folder, genome, fdr_thr, numberclust)
		sink()
		setwd(current_path)
	}

	# Inter IGV html file generation
	if (file.exists(paste0(diffdir,"/fdr_result/differential.inter_sample_group.",pcgrp,".bedGraph"))) {
		compartment_file <- read.table(paste0(diffdir,"/fdr_result/differential.inter_sample_group.",pcgrp,".bedGraph"), h=T, as.is=T)
		compartment_file <- compartment_file[order(compartment_file$chr, compartment_file$start),c("chr","start","end",prefix_master,"replicate_wt","sample_maha","padj","dist_clust")]
		print (head(compartment_file))
		compartment_file$padj[compartment_file$padj == 0] <- min(compartment_file$padj[compartment_file$padj > 0])
		compartment_file$padj <- -log10(compartment_file$padj)
		write.table(compartment_file[,1:(ncol(compartment_file))], file=paste0(diffdir,"/viz/files/inter_compartment.bedGraph"), row.names=F, col.names=T, sep="\t", quote=F)

		inter_pcori <- list.files(paste0(diffdir,"/",pcgrp), pattern="inter_sample", full.names=T)
		pcori_files <- list()
		for(i in 1:length(inter_pcori)) {
			pcori_files[[i]] <- read.table(inter_pcori[i], h=T, as.is=T)
		}

		pcori_files <- do.call(rbind, pcori_files)
		pcori_files <- pcori_files[order(pcori_files$chr, pcori_files$start),]
		for(i in 1:(ncol(pcori_files) - 3)) {
			write.table(pcori_files[,c(1:3,(i+3))], file=paste0(diffdir,"/viz/files/inter_",colnames(pcori_files)[i+3],"_PC.bedGraph"), row.names=F, col.names=F, sep="\t", quote=F)
		}

		if (file.exists(paste0(diffdir,"/fdr_result/inter_sample_group.subcompartments.bedGraph"))) {
			subcomp <- read.table(paste0(diffdir,"/fdr_result/inter_sample_group.subcompartments.bedGraph"), h=T, as.is=T)
			for(j in 1:length(prefix_master)) {
				subcomp_temp <- subcomp[,c("chr","start","end",paste0(prefix_master[j],".score"))]
				colnames(subcomp_temp)  <- c("Chromosome","Start","End","Segment_Mean")
				subcomp_temp[,"Sample"] <- prefix_master[j]
				subcomp_temp[,"Num_Probes"] <- 1
				subcomp_temp <- subcomp_temp[,c("Sample","Chromosome","Start","End","Num_Probes","Segment_Mean")]
				write.table(subcomp_temp, file=paste0(diffdir,"/viz/files/inter_",prefix_master[j],"_subcomp.seg"), row.names=F, col.names=T, sep="\t", quote=F)
			}
		}

		file <- list()
		file[[1]] <- data.frame(
					file=paste0("./files/inter_compartment.bedGraph"),
					name="Selected",
					group="compartment",
					color="234,100,0:122,16,180"
				)

		comp_colA <- t(col2rgb(colorspace::diverge_hcl(length(prefix_master), c=100, l=c(50,90), power=0.5)))
		comp_colB <- t(col2rgb(colorspace::heat_hcl(length(prefix_master), h=c(80,-100), l=c(75,40), c=c(40,80), power=0.5)))
		#k <- length(prefix_master)
		#l <- 2
		#for(i in 1:length(prefix_master)) {
		#	data_rep <- data[data$prefix.master== as.character(prefix_master[i]),] 
		#	print (data_rep)
		#	if (comp_colA[k,1] > comp_colB[k,1]) {
		#		red_col <- paste0(comp_colA[k,], collapse=",")
		#		blu_col <- paste0(comp_colB[k,], collapse=",")
		#	} else {
		#		red_col <- paste0(comp_colB[k,], collapse=",")
		#		blu_col <- paste0(comp_colA[k,], collapse=",")
		#	}
		#	for(j in 1:nrow(data_rep)) {
		#		file[[l]] <- data.frame(
		#					file=paste0("./files/inter_",data_rep$prefix[j],"_PC.bedGraph"),
		#					name=data_rep$prefix[j],
		#					group="bedGraph",
		#					color=paste0(red_col,":",blu_col)
		#				)
		#		l <- l + 1
		#	}
		#	k <- k - 1
		#}
		
		if (file.exists(paste0(diffdir,"/fdr_result/inter_sample_group.subcompartments.bedGraph"))) {
			for(j in 1:length(prefix_master)) {
				file[[length(file)+j]] <- data.frame(
							file=paste0("./files/inter_",prefix_master[j],"_subcomp.seg"),
							name=paste0(prefix_master[j],".subcomp"),
							group="seg",
							color="0:0"
						)
			}
		}
		file <- as.data.frame(do.call(rbind, file), stringsAsFactors=F)
		print (file)

		current_path <- getwd()
		setwd(paste0(diffdir,"/viz/"))
		folder <- "vizIGV_inter"
		if (!file.exists(folder)) {
 			dir.create(folder)
		}
		sink(paste0(folder,"/inter_igv.html"))
		htmlheader()
		htmlbody(file, folder, genome, fdr_thr, numberclust)
		sink()
		setwd(current_path)
	}
}

### END ###

########################### Gene enrichment ###########################

geneEnrichment <- function(data, diffdir, genome, exclA=T, region="anchor", pcgrp="pcQnm", interaction="intra", pcscore=F, compare=F, cells="X,Y", fdr_thr=0.05) {

	print (prefix_master)
	prefix_master <- unique(as.character(data$prefix.master))
	if (compare) {
		prefix_master <- strsplit(cells,",")[[1]]
	}

	folder <- paste0("DifferentialResult/",diffdir,"/geneEnrichment")
	if(!dir.exists(folder)) {
		dir.create(folder)
	}
	if (compare == TRUE) {
		folder <- paste0("DifferentialResult/",diffdir,"/geneEnrichment/comparison_",paste0(strsplit(cells,",")[[1]],collapse="_vs_"))
		if (!dir.exists(folder)) {
			dir.create(paste0("DifferentialResult/",diffdir,"/geneEnrichment/comparison_",paste0(strsplit(cells,",")[[1]],collapse="_vs_")))
		}
	}

	compartment_score <- read.table(paste0("DifferentialResult/",diffdir,"/fdr_result/differential.",interaction,"_sample_group.",pcgrp,".bedGraph"), h=T, as.is=T)
	diff_compartments <- read.table(paste0("DifferentialResult/",diffdir,"/fdr_result/differential.",interaction,"_sample_group.Filtered.",pcgrp,".bedGraph"), h=T, as.is=T)
	rownames(compartment_score) <- paste0(compartment_score$chr,"_",compartment_score$start)
	rownames(diff_compartments) <- paste0(diff_compartments$chr,"_",diff_compartments$start)
	resolution <- compartment_score$end[1] - compartment_score$start[1]
	resolution <- as.integer(resolution)

	print (paste0(genome,"_",resolution,"_goldenpathData"))
	if (!dir.exists(paste0(genome,"_",resolution,"_goldenpathData"))) {
		cat(paste0(genome,"_",resolution,"_goldenpathData does not exist, creating!"))

		## Data download and processing 
    	folder_genome <- paste0(genome,"_",as.integer(resolution),"_goldenpathData")
    	dir.create(folder_genome)
    	folder_genome <- normalizePath(folder_genome)

  		current_path <- getwd()
		setwd(folder_genome)
		datadownload(genome)
  		setwd(current_path)

 		if (!file.exists(paste0(folder_genome,"/",genome,".fa"))) {
   			cat ("Unzipping ",paste0(folder_genome,"/",genome,".fa.gz"),"\n")
   			system(paste0("gunzip -c ",folder_genome,"/",genome,".fa.gz > ",folder_genome,"/",genome,".fa"), wait=T)
		}
		if (!file.exists(paste0(folder_genome,"/",genome,".tss.bed"))) {
   			cmd <- paste0("gunzip -c ",folder_genome,"/",genome,".refGene.gtf.gz |awk -v OFS='\\t' '{if($3==\"transcript\"){if($7==\"+\"){print $1,$4,$4+1}else{print $1,$5-1,$5}}}' |grep -v \"alt\" |grep -v \"random\" |sort |uniq |sort -k 1,1 -k2,2n > ",folder_genome,"/",genome,".tss.bed")
   			cat ("Running ",cmd,"\n")
   			system(cmd, wait=T)
		}
		if (!file.exists(paste0(folder_genome,"/",genome,".binned.bed"))) {
  			cmd <- paste0(.findExecutable("bedtools")," makewindows -g ",folder_genome,"/",genome,".chrom.sizes -w ",as.integer(resolution)," > ",folder_genome,"/",genome,".binned.bed")
   			cat ("Running ",cmd,"\n")
   			system(cmd, wait=T)
		}
		if (!file.exists(paste0(folder_genome,"/",genome,".GCpt.bedGraph"))) {
  			cmd <- paste0(.findExecutable("bedtools")," nuc -fi ",folder_genome,"/",genome,".fa -bed ",folder_genome,"/",genome,".binned.bed |grep -v \"#\" |awk -v OFS='\\t' '{print $1,$2,$3,$5}' |grep -v \"alt\" |grep -v \"random\" |sort -k 1,1 -k2,2n > ",folder_genome,"/",genome,".GCpt.bedGraph")
  			cat ("Running ",cmd,"\n")
   			system(cmd, wait=T)
		}
		if (!file.exists(paste0(folder_genome,"/",genome,".GCpt.tss.bedGraph"))) {
			cmd <- paste0(.findExecutable("bedtools")," map -a ",folder_genome,"/",genome,".GCpt.bedGraph -b ",folder_genome,"/",genome,".tss.bed -c 1 -o count -null 0 > ",folder_genome,"/",genome,".GCpt.tss.bedGraph")
    		cat ("Running ",cmd,"\n")
    		system(cmd, wait=T)
  		}

	}
	refgene <- read.table(paste0(genome,"_",resolution,"_goldenpathData/",genome,".refGene.gtf.gz"), h=F, sep="\t", fill=NA)
	refgene_bed <- refgene[refgene$V3=="transcript",c(1,4,5)]
	colnames(refgene_bed) <- c("chr","start","end")
	refgene_bed[,"gene"]  <- gsub(";","",do.call(rbind,strsplit(as.character(refgene[refgene$V3=="transcript",]$V9)," "))[,2])
	refgene_bed <- unique(refgene_bed)
	refgene_bed <- refgene_bed[order(refgene_bed$chr, refgene_bed$start),]
	write.table(refgene_bed, file=paste0(folder,"/refgene.bed"), row.names=F, col.names=F, sep="\t", quote=F)

	if (region == "anchor" | region == "interactor" | region == "both") {
		if (exclA == FALSE & pcgrp == "pcOri") {
			compartment_score_pcQnm <- read.table(paste0("DifferentialResult/",diffdir,"/fdr_result/differential.",interaction,"_sample_group.pcQnm.bedGraph"), h=T, as.is=T)
			diff_compartments_pcQnm <- read.table(paste0("DifferentialResult/",diffdir,"/fdr_result/differential.",interaction,"_sample_group.Filtered.pcQnm.bedGraph"), h=T, as.is=T)
			rownames(compartment_score_pcQnm) <- paste0(compartment_score_pcQnm$chr,"_",compartment_score_pcQnm$start)
			rownames(diff_compartments_pcQnm) <- paste0(diff_compartments_pcQnm$chr,"_",diff_compartments_pcQnm$start)
			for(i in 1:length(prefix_master)) {
				if(!dir.exists(paste0(folder,"/",prefix_master[i],"_geneEnrichment"))) {
					dir.create(paste0(folder,"/",prefix_master[i],"_geneEnrichment"))
				}
				diff_compartments_sample_A <- diff_compartments[,c("chr","start","end",prefix_master[i])]
				colnames(diff_compartments_sample_A)[4] <- "sample"
				diff_compartments_sample_A <- diff_compartments_sample_A[diff_compartments_sample_A$sample > 0,]
				if (pcscore) {
					diff_compartments_pcQnm_sample <- diff_compartments_pcQnm[,c("chr","start","end",prefix_master[i])]
					colnames(diff_compartments_pcQnm_sample)[4] <- "sample"
					if (length(prefix_master) > 2) {
						diff_compartments_pcQnm_rest <- diff_compartments_pcQnm[,c(prefix_master[which(prefix_master != prefix_master[i])])]
						diff_compartments_pcQnm_sample[,"sample_rest"] <- apply(diff_compartments_pcQnm_rest, 1, max)
					} else {
						diff_compartments_pcQnm_sample[,"sample_rest"] <- diff_compartments_pcQnm[,c(prefix_master[which(prefix_master != prefix_master[i])])]
					}
					diff_compartments_pcQnm_sample <- diff_compartments_pcQnm_sample[diff_compartments_pcQnm_sample$sample > diff_compartments_pcQnm_sample$sample_rest,]
					diff_compartments_sample_A <- na.omit(diff_compartments_sample_A[rownames(diff_compartments_pcQnm_sample),])
					write.table(diff_compartments_sample_A[,1:3], file=paste0(folder,"/",prefix_master[i],"_geneEnrichment/",prefix_master[i],"_Diff_A_compartments.bedGraph"), row.names=F, col.names=F, sep="\t", quote=FALSE)
				} else {
					write.table(diff_compartments_sample_A[,1:3], file=paste0(folder,"/",prefix_master[i],"_geneEnrichment/",prefix_master[i],"_Diff_A_compartments.bedGraph"), row.names=F, col.names=F, sep="\t", quote=FALSE)
				}
			}
		} else if (exclA == TRUE & (pcgrp == "pcOri" | pcgrp == "pcQnm")) {
			for(i in 1:length(prefix_master)) {
				if(!dir.exists(paste0(folder,"/",prefix_master[i],"_geneEnrichment"))) {
					dir.create(paste0(folder,"/",prefix_master[i],"_geneEnrichment"))
				}
				diff_compartments_sample_A <- diff_compartments[,c("chr","start","end",prefix_master[i])]
				colnames(diff_compartments_sample_A)[4] <- "sample"
				if (length(prefix_master) > 2) {
					diff_compartments_sample_rest <- diff_compartments[,c(prefix_master[which(prefix_master != prefix_master[i])])]
					diff_compartments_sample_A[,"sample_rest"] <- apply(diff_compartments_sample_rest, 1, max)
				} else {
					diff_compartments_sample_A[,"sample_rest"] <- diff_compartments[,c(prefix_master[which(prefix_master != prefix_master[i])])]
				}
				if (pcscore) {
					cat ("When exclA option is set as TRUE, pcscore option will be treated as FALSE\n")
				}
				diff_compartments_sample_A <- diff_compartments_sample_A[diff_compartments_sample_A$sample > 0 & diff_compartments_sample_A$sample_rest < 0,]
				write.table(diff_compartments_sample_A[,1:3], file=paste0(folder,"/",prefix_master[i],"_geneEnrichment/",prefix_master[i],"_Diff_A_compartments.bedGraph"), row.names=F, col.names=F, sep="\t", quote=FALSE)
			}
		}

		if (region == "anchor") {
			for(i in 1:length(prefix_master)) {
				cmd <- paste0("bedtools intersect -a ",folder,"/",prefix_master[i],"_geneEnrichment/",prefix_master[i],"_Diff_A_compartments.bedGraph -b ",folder,"/refgene.bed -wb |awk '{print $NF}' |sort |uniq")
				gene_list <- system(cmd, wait=T, intern=T)
				gene_list <- paste0("\"",gene_list,"\"",collapse=",")
				cat(paste0("#!/bin/sh\ncurl -H 'Content-Type: text/json' -d \'{\"Symbols\":[",gene_list,"]}' https://toppgene.cchmc.org/API/lookup > ",folder,"/",prefix_master[i],"_geneEnrichment/lookup_result.json"), file=paste0(folder,"/",prefix_master[i],"_geneEnrichment/lookup.sh"))
				system(paste0("chmod 755 ",folder,"/",prefix_master[i],"_geneEnrichment/lookup.sh"), wait=T)
				system(paste0("./",folder,"/",prefix_master[i],"_geneEnrichment/lookup.sh"), wait=T)

				print (paste0(folder,"/",prefix_master[i],"_geneEnrichment/lookup_result.json"))
				converted_genes <- rjson::fromJSON(file=paste0(folder,"/",prefix_master[i],"_geneEnrichment/lookup_result.json"))
				gene_list_tmp <- unique(as.character(unlist(converted_genes)[which(names(unlist(converted_genes)) == "Genes.Entrez")]))
				gene_list <- paste0(gene_list_tmp,collapse=",")
				cat(paste0("#!/bin/sh\ncurl -H 'Content-Type: text/json' -d \'{\"Genes\":[",gene_list,"]}' https://toppgene.cchmc.org/API/enrich > ",folder,"/",prefix_master[i],"_geneEnrichment/gene_enrich.json"), file=paste0(folder,"/",prefix_master[i],"_geneEnrichment/enrich.sh"))
				system(paste0("chmod 755 ",folder,"/",prefix_master[i],"_geneEnrichment/enrich.sh"), wait=T)
				system(paste0("./",folder,"/",prefix_master[i],"_geneEnrichment/enrich.sh"), wait=T)

				if (file.size(paste0(folder,"/",prefix_master[i],"_geneEnrichment/gene_enrich.json")) > 110) {
					enrich_list   <- list()
					enrich_result <- rjson::fromJSON(file = paste0(folder,"/",prefix_master[i],"_geneEnrichment/gene_enrich.json"))
					terms_enrched <- length(do.call(c,enrich_result))
					for(j in 1:terms_enrched) {
						enrich_list[[j]] <- as.data.frame(do.call(cbind,do.call(c,enrich_result)[[j]]))[1,1:11]
						enrich_list[[j]][,"Genes"] <- paste0(unlist(do.call(cbind,do.call(c,enrich_result)[[j]]$Genes)[2,]), collapse=",")
					}
					enrich_list <- do.call(rbind, enrich_list)
					enrich_list <- as.data.frame(enrich_list)
					enrich_list <- enrich_list[enrich_list$QValueFDRBH < 0.05,]
					write.table(apply(enrich_list,2,as.character), file=paste0(folder,"/",prefix_master[i],"_geneEnrichment/",prefix_master[i],"_geneEnrichment.",region,".txt"), row.names=F, col.names=T, sep="\t", quote=F)
					write.table(gene_list_tmp, file=paste0(folder,"/",prefix_master[i],"_geneEnrichment/",prefix_master[i],"_geneList.",region,".txt"), row.names=F, col.names=F, sep="\t", quote=F)
				} else {
					cat("There are way to many genes overlapping with the A compartment, the command line option has return an error!\nPlease mannually upload the gene list.")
					write.table(gene_list_tmp, file=paste0(folder,"/",prefix_master[i],"_geneEnrichment/",prefix_master[i],"_geneList.",region,".txt"), row.names=F, col.names=F, sep="\t", quote=F)
				}
			}
		}
	} 

	if (region == "interactor" | region == "both") {
		if (interaction == "inter") {
			stop ("Function enrichment of interactor region is not possible for inter compartments differential loops\n")
		} else {
			bedpe <- read.table(paste0("DifferentialResult/",diffdir,"/fdr_result/differential.intra_compartmentLoops.bedpe"), h=T, as.is=T)
			bedpe <- bedpe[order(bedpe$chr1, bedpe$start1),]
			bedpe[,"id1"] <- paste0(bedpe$chr1,"_",bedpe$start1)
			bedpe[,"id2"] <- paste0(bedpe$chr2,"_",bedpe$start2)
		
			for(i in 1:length(prefix_master)) {
				compartment_score_sample <- compartment_score[,c("chr","start","end",prefix_master[i],"padj")]
				diff_compartments_sample <- read.table(paste0(folder,"/",prefix_master[i],"_geneEnrichment/",prefix_master[i],"_Diff_A_compartments.bedGraph"), h=F, as.is=T)
				colnames(diff_compartments_sample) <- c("chr","start","end")
				rownames(diff_compartments_sample) <- paste0(diff_compartments_sample$chr,"_",diff_compartments_sample$start)
				bedpe[,c("pcscore1","padj1")] <- na.omit(compartment_score_sample[as.character(bedpe$id1),c(4,5)])
				bedpe[,c("pcscore2","padj2")] <- na.omit(compartment_score_sample[as.character(bedpe$id2),c(4,5)])
				bedpe[,"class1"] <- "no"
				bedpe[,"class2"] <- "no"
				bedpe[bedpe$padj1 < fdr_thr, "class1"] <- "yes"
				bedpe[bedpe$padj2 < fdr_thr, "class2"] <- "yes"
				bedpe <- bedpe[(bedpe$class1 == "yes" & bedpe$class2 == "no") | (bedpe$class1 == "no" & bedpe$class2 == "yes"),]
				k <- 0
				bedpe_filt <- list()
				for(j in 1:nrow(bedpe)) {
					if (bedpe$class1[j] == "yes") {
						id <- bedpe$id1[j]
					} else if (bedpe$class2[j] == "yes") {
						id <- bedpe$id2[j]
					}
					if (length(which(id == rownames(diff_compartments_sample))) == 1) {
						k <- k + 1
						bedpe_filt[[k]] <- bedpe[j,]
					}
				}
				bedpe_filt <- do.call(rbind, bedpe_filt)
				bed_collapse_A <- data.frame(chr=c(bedpe_filt$chr1,bedpe_filt$chr2), start=c(bedpe_filt$start1,bedpe_filt$start2), end=c(bedpe_filt$end1,bedpe_filt$end2), pcscore=c(bedpe_filt$pcscore1,bedpe_filt$pcscore2), class=c(bedpe_filt$class1,bedpe_filt$class2)) 
				bed_collapse_A <- bed_collapse_A[bed_collapse_A$class == "no",]
				bed_collapse_A <- unique(bed_collapse_A[bed_collapse_A$pcscore > 0,])
				write.table(bed_collapse_A[,1:3], file=paste0(folder,"/",prefix_master[i],"_geneEnrichment/",prefix_master[i],"_Diff_A_compartments_interactor.bedGraph"), row.names=F, col.names=F, sep="\t", quote=FALSE)
			}

			if (region == "interactor") {
				for(i in 1:length(prefix_master)) {
					cmd <- paste0("bedtools intersect -a ",folder,"/",prefix_master[i],"_geneEnrichment/",prefix_master[i],"_Diff_A_compartments_interactor.bedGraph -b ",folder,"/refgene.bed -wb |awk '{print $NF}' |sort |uniq")
					gene_list <- system(cmd, wait=T, intern=T)
					gene_list <- paste0("\"",gene_list,"\"",collapse=",")
					cat(paste0("#!/bin/sh\ncurl -H 'Content-Type: text/json' -d \'{\"Symbols\":[",gene_list,"]}' https://toppgene.cchmc.org/API/lookup > ",folder,"/",prefix_master[i],"_geneEnrichment/lookup_result.json"), file=paste0(folder,"/",prefix_master[i],"_geneEnrichment/lookup.sh"))
					system(paste0("chmod 755 DifferentialResult/",diffdir,"/geneEnrichment/",prefix_master[i],"_geneEnrichment/lookup.sh"), wait=T)
					system(paste0("./",folder,"/",prefix_master[i],"_geneEnrichment/lookup.sh"), wait=T)

					print (paste0(folder,"/",prefix_master[i],"_geneEnrichment/lookup_result.json"))
					converted_genes <- rjson::fromJSON(file=paste0(folder,"/",prefix_master[i],"_geneEnrichment/lookup_result.json"))
					gene_list_tmp <- unique(as.character(unlist(converted_genes)[which(names(unlist(converted_genes)) == "Genes.Entrez")]))
					gene_list <- paste0(gene_list_tmp,collapse=",")
					cat(paste0("#!/bin/sh\ncurl -H 'Content-Type: text/json' -d \'{\"Genes\":[",gene_list,"]}' https://toppgene.cchmc.org/API/enrich > ",folder,"/",prefix_master[i],"_geneEnrichment/gene_enrich.json"), file=paste0(folder,"/",prefix_master[i],"_geneEnrichment/enrich.sh"))
					system(paste0("chmod 755 ",folder,"/",prefix_master[i],"_geneEnrichment/enrich.sh"), wait=T)
					system(paste0("./",folder,"/",prefix_master[i],"_geneEnrichment/enrich.sh"), wait=T)

					if (file.size(paste0(folder,"/",prefix_master[i],"_geneEnrichment/gene_enrich.json")) > 110) {
						enrich_list   <- list()
						enrich_result <- rjson::fromJSON(file = paste0(folder,"/",prefix_master[i],"_geneEnrichment/gene_enrich.json"))
						terms_enrched <- length(do.call(c,enrich_result))
						for(j in 1:terms_enrched) {
							enrich_list[[j]] <- as.data.frame(do.call(cbind,do.call(c,enrich_result)[[j]]))[1,1:11]
							enrich_list[[j]][,"Genes"] <- paste0(unlist(do.call(cbind,do.call(c,enrich_result)[[j]]$Genes)[2,]), collapse=",")
						}
						enrich_list <- do.call(rbind, enrich_list)
						enrich_list <- as.data.frame(enrich_list)
						enrich_list <- enrich_list[enrich_list$QValueFDRBH < 0.05,]
						write.table(apply(enrich_list,2,as.character), file=paste0(folder,"/",prefix_master[i],"_geneEnrichment/",prefix_master[i],"_geneEnrichment.",region,".txt"), row.names=F, col.names=T, sep="\t", quote=F)
						write.table(gene_list_tmp, file=paste0(folder,"/",prefix_master[i],"_geneEnrichment/",prefix_master[i],"_geneList.",region,".txt"), row.names=F, col.names=F, sep="\t", quote=F)
					} else {
						cat("There are way to many genes overlapping with the A compartment, the command line option has return an error!\nPlease mannually upload the gene list.")
						write.table(gene_list_tmp, file=paste0(folder,"/",prefix_master[i],"_geneEnrichment/",prefix_master[i],"_geneList.",region,".txt"), row.names=F, col.names=F, sep="\t", quote=F)
					}
				}
			}

			if (region == "both") {
				for(i in 1:length(prefix_master)) {
					diff_compartment_sample <- read.table(paste0(folder,"/",prefix_master[i],"_geneEnrichment/",prefix_master[i],"_Diff_A_compartments.bedGraph"), h=F, as.is=T)
					diff_interactors_sample <- read.table(paste0(folder,"/",prefix_master[i],"_geneEnrichment/",prefix_master[i],"_Diff_A_compartments_interactor.bedGraph"), h=F, as.is=T)
					colnames(diff_compartment_sample) <- c("chr","start","end")
					colnames(diff_interactors_sample) <- c("chr","start","end")
					diff_both <- as.data.frame(rbind(diff_compartment_sample, diff_interactors_sample))
					diff_both <- diff_both[order(diff_both$chr, diff_both$start),]
					write.table(diff_both, file=paste0(folder,"/",prefix_master[i],"_geneEnrichment/",prefix_master[i],"_Diff_A_both.bedGraph"), row.names=F, col.names=F, sep="\t", quote=FALSE)

					cmd <- paste0("bedtools intersect -a ",folder,"/",prefix_master[i],"_geneEnrichment/",prefix_master[i],"_Diff_A_both.bedGraph -b ",folder,"/refgene.bed -wb |awk '{print $NF}' |sort |uniq")
					gene_list <- system(cmd, wait=T, intern=T)
					gene_list <- paste0("\"",gene_list,"\"",collapse=",")
					cat(paste0("#!/bin/sh\ncurl -H 'Content-Type: text/json' -d \'{\"Symbols\":[",gene_list,"]}' https://toppgene.cchmc.org/API/lookup > ",folder,"/",prefix_master[i],"_geneEnrichment/lookup_result.json"), file=paste0(folder,"/",prefix_master[i],"_geneEnrichment/lookup.sh"))
					system(paste0("chmod 755 ",folder,"/",prefix_master[i],"_geneEnrichment/lookup.sh"), wait=T)
					system(paste0("./",folder,"/",prefix_master[i],"_geneEnrichment/lookup.sh"), wait=T)

					print (paste0(folder,"/",prefix_master[i],"_geneEnrichment/lookup_result.json"))
					converted_genes <- rjson::fromJSON(file=paste0(folder,"/",prefix_master[i],"_geneEnrichment/lookup_result.json"))
					gene_list_tmp <- unique(as.character(unlist(converted_genes)[which(names(unlist(converted_genes)) == "Genes.Entrez")]))
					gene_list <- paste0(gene_list_tmp,collapse=",")
					cat(paste0("#!/bin/sh\ncurl -H 'Content-Type: text/json' -d \'{\"Genes\":[",gene_list,"]}' https://toppgene.cchmc.org/API/enrich > ",folder,"/",prefix_master[i],"_geneEnrichment/gene_enrich.json"), file=paste0(folder,"/",prefix_master[i],"_geneEnrichment/enrich.sh"))
					system(paste0("chmod 755 ",folder,"/",prefix_master[i],"_geneEnrichment/enrich.sh"), wait=T)
					system(paste0("./",folder,"/",prefix_master[i],"_geneEnrichment/enrich.sh"), wait=T)

					if (file.size(paste0(folder,"/",prefix_master[i],"_geneEnrichment/gene_enrich.json")) > 110) {
						enrich_list   <- list()
						enrich_result <- rjson::fromJSON(file = paste0(folder,"/",prefix_master[i],"_geneEnrichment/gene_enrich.json"))
						terms_enrched <- length(do.call(c,enrich_result))
						for(j in 1:terms_enrched) {
							enrich_list[[j]] <- as.data.frame(do.call(cbind,do.call(c,enrich_result)[[j]]))[1,1:11]
							enrich_list[[j]][,"Genes"] <- paste0(unlist(do.call(cbind,do.call(c,enrich_result)[[j]]$Genes)[2,]), collapse=",")
						}
						enrich_list <- do.call(rbind, enrich_list)
						enrich_list <- as.data.frame(enrich_list)
						enrich_list <- enrich_list[enrich_list$QValueFDRBH < 0.05,]
						write.table(apply(enrich_list,2,as.character), file=paste0(folder,"/",prefix_master[i],"_geneEnrichment/",prefix_master[i],"_geneEnrichment.",region,".txt"), row.names=F, col.names=T, sep="\t", quote=F)
						write.table(gene_list_tmp, file=paste0(folder,"/",prefix_master[i],"_geneEnrichment/",prefix_master[i],"_geneList.",region,".txt"), row.names=F, col.names=F, sep="\t", quote=F)
					} else {
						cat("There are way to many genes overlapping with the A compartment, the command line option has return an error!\nPlease mannually upload the gene list.")
						write.table(gene_list_tmp, file=paste0(folder,"/",prefix_master[i],"_geneEnrichment/",prefix_master[i],"_geneList.",region,".txt"), row.names=F, col.names=F, sep="\t", quote=F)
					}
				}
			}
		}
	}
}

### END ###

#Get the options 
option_list = list(
	make_option(c("--file"), type="character", help="an input.txt file with the following columns, e.g.\n
	\t\t<mat>\t<bed>\t<prefix_replicate>\t<prefix>
	\t\tmat.txt\tmat.bed\tname_R1\tname\n
	Note: Do not use '-'/dash and '.'/dot in the prefix_replicate and prefix name\n 
	Here mat.txt file is an upper triangular Hi-C sparse matrix
	It should have the following columns\n
	\t\t<indexA> <indexB> <count>\n
 	\t\t1 1 300
 	\t\t1 2 30
 	\t\t1 3 10
 	\t\t2 2 200
 	\t\t2 3 20
 	\t\t3 3 200
 	\t\t....\n
 	mat.bed is bed file with index information
 	It should have the following columns\n
 	\t\t<chr>\t<start>\t<end>\t<index>\n
 	\t\tchr1\t0\t40000\t1
 	\t\tchr1\t40000\t80000\t2
 	\t\tchr1\t80000\t120000\t3
 	\t\t....\n
 	Optionally it can have an extra column that will denote if the bin is a blacklisted (1) or not (0)\n
 	\t\t<chr>\t<start>\t<end>\t<index>\t<blacklisted>\n
 	\t\tchr1\t0\t40000\t1\t0
 	\t\tchr1\t40000\t80000\t2\t1
 	\t\tchr1\t80000\t120000\t3\t1
 	\t\t....\n
 	dcHiC will omit the blacklisted regions during PCA calculation\n"),

 	make_option(c("--pcatype"), type="character", default="cis", help="pcatype will allow a user to perform different types of step wise analyzis on the Hi-C data. 
 		Available options: cis/trans/select/analyze/subcomp/dloop/viz, for more details check the small description and example commands to run the analysis\n 

 		cis     : Perform PCA on cis interaction matrix, e.g.
		Rscript ./dchicf.r --file input.txt --pcatype cis --dirovwt T

 		trans   : Perform PCA on trans interaction matrix, e.g.
 		Rscript ./dchicf.r --file input.txt --pcatype trans --dirovwt T

 		select  : Selection of best PC for downstream analysis [Must after cis or trans step], e.g.
 		Rscript ./dchicf.r --file input.txt --pcatype select --dirovwt T --genome mm10 

 		analyze : Perform differential PCA analysis on the selected PC's across the samples [Must after select step], e.g. 
 		Rscript ./dchicf.r --file input.txt --pcatype analyze --dirovwt T --diffdir conditionA_vs_conditionB

 		subcomp : Assigning sub-compartments based in the PC magnitude values using HMM segmentation [Optional step to classify the compartments based on their magnitude], e.g. 
 		Rscript ./dchicf.r --file input.txt --pcatype subcomp --diffdir conditionA_vs_conditionB

 		fithic  : Before calling the differential loop (dloop) it is advisable to run the fithic part [Must before the differential loop calling part], e.g. 
 		Rscript ./dchicf.r --file input.txt --pcatype fithic --dirovwt T --diffdir conditionA_vs_conditionB --maxd 10e6 --fithicpath '/path/to/fithic.py' --pythonpath '/path/to/python3.x/python'


 		dloop   : Find differential loops anchored in at least one of the differential compartments across the samples [Optional but recommended step to find differential loops originating from differential compartments], e.g. 
 		Rscript ./dchicf.r --file input.txt --pcatype dloop --diffdir conditionA_vs_conditionB --maxd 10e6 --minc 0 

 		viz     : Generate IGV vizualization html file [Must inorder to vizualize the result]
 		Rscript ./dchicf.r --file input.txt --pcatype viz --diffdir conditionA_vs_conditionB --genome mm10 --pcgroup pcQnm

 		enrich  : Perform gene enrichment analysis [Gene set enrichment of differential genes/loops overlapping genes], e.g.
 		Rscript ./dchicf.r --file input.txt --pcatype enrich --genome mm10 --diffdir conditionA_vs_conditionB --exclA F --region both --pcgroup pcQnm --interaction intra --pcscore F --compare F
 		\n"),

 	make_option(c("--ebackgrnd"), type="integer", default=2, help="Expected genomic background to calculate Observed/Expected ratio
 		1 : Scheme 1 will use the full genomic background to calculate the expected counts 
 		2 : Scheme 2 will chromosome wise background to estimate the expected counts [default]\n"),
 	
 	make_option(c("--pc"), type="integer", default=2, help="Number of PC's to be written 
 		[default 2]\n"),
 	
 	make_option(c("--rowmrge"), type="integer", default=10000, help="This is the number of consecutive rows to be merged in the trans matrix file and add the column count values to reduce the overall size 
 		[rows=trans bins, cols=cis bins, default=10,000]\n"),
 	
 	make_option(c("--genome"), type="character", default=NA, help="Genome name, e.g. hg38/hg19/mm10/mm9/[custom]
 		[default NONE, should be provided by the user]\n"),

	make_option(c("--gfolder"), type="character", default=NA, help="Genome folder path, e.g. <genome>_<resolution>_goldenpathData folder (This folder should contain three files - <genome>.fa, <genome>.tss.bed and <genome>.chrom.sizes files) 
                [default NA]\n"),

 	make_option(c("--bdgfile"), type="character", default=NA, help="Additional bedGraph files to be visualized in the html file. Optional for viz step. This is a single three column txt file\n
  		<path to bedGraph>\t<name>\t<type>\t<pos color:neg color>\n
  		The pos and neg color is the RGB code will represent the respective postive and negative score distribution in the bedGraph file
  		If you have only postive score, keep the neg part as 0, e.g. 250,0,0:0\n"),
 	
 	make_option(c("--rzsc"), type="numeric", default=2, help="Replicate zscore threshold. Increasing the value makes dcHiC non-specific towards replicate-wise PC variation  
 		[default 2]\n"),

 	make_option(c("--szsc"), type="numeric", default=0, help="Sample zscore threshold. Increasing the value makes dcHiC specific towards finding large PC difference across sample
 		[default 0]\n"),

 	make_option(c("--fdr"), type="numeric", default=0.05, help="Report only above this FDR threshold
 		[default 0.05]\n"),

 	make_option(c("--distclust"), type="numeric", default=-1, help="Distance threshold for clustering close differential regions. This will help to remove lone differential compartments in the chromosome. A distclust value of 0 will cluster differential compartments that are adjacent to each other. Higher the value further apart the differential compartments are
 		[default -1 i.e. no clustering of result]\n"),

 	make_option(c("--numberclust"), type="numeric", default=1, help="Minimum number of clustered regions within the predefined distclust threshold (when it is set above default i.e. -1)
 		[default 1]\n"),

 	make_option(c("--subnum"), type="integer", default=6, help="Total number of sub-compartments to find 
 		[default 6]\n"),

 	 make_option(c("--minc"), type="integer", default=5, help="Minimum interaction count to be retained during differential loop calling 
 		[default 6]\n"),

 	 make_option(c("--maxd"), type="integer", default=2e6, help="Maximum distance (bp) between interactions to be kept during differential loop calling 
 		[default 2e6]\n"),

	make_option(c("--exclA"), type="logical", default=TRUE, help="Setting this option to TRUE will allow dcHiC to search genes that are overlapping with differential-A regions in a sample and B in all others [default TRUE]\n"),

	make_option(c("--pcscore"), type="logical", default=FALSE, help="Setting this option to TRUE will allow dcHiC to search genes that are overlapping with differential-A regions with higher pc score than rest of the samples [default FALSE]\n"),

	make_option(c("--interaction"), type="character", default="intra", help="dcHiC can perform functional enrichment over both cis/intra and trans/inter compartments. This will allow the user to choose the compartment type [default intra]\n"),

	make_option(c("--compare"), type="logical", default=F, help="dcHiC can perform functional enrichment between a pair of defined sample. Set this option to TRUE when user needs to compare function enrichment of A-compartment genes bwteen a pair instead all the samples [default FALSE]\n"),

	make_option(c("--cells"), type="character", default="X,Y", help="When the compare argument is set to TRUE then dcHiC expects the name <prefix> of sample pairs provided in the input.txt file [example prefix1,prefix2]\n"),

	make_option(c("--region"), type="character", default="both", help="dcHiC finds differential compartments (anchors) and then differential interactions (interactors) anchored in differential compartments
		Setting this option to 'anchor' will only search genes overlapping with differential compartments. When set to 'interactors', the genes overlaping with the ineracting regions will be scanned
		[anchor/interactor/both default both]\n"),

	make_option(c("--pcgroup"), type="character", default="pcQnm", help="Setting the option to pcQnm will direct dcHiC to use quantile normalized PC values to create IGV html file. When set to pcOri, dcHiC will use the original PC values in IGV html file (viz step)
		Gene enrichment step will also be effected by the option accordingly.
		[pcQnm/pcOri default pcQnm]\n"),

 	make_option(c("--diffdir"), type="character", default="sample", help="Folder name under which all the differential results will be saved
 		[default sample]\n"),

 	make_option(c("--dirovwt"), type="logical", default=FALSE, help="Overwrite the existing directories and files 
 		[default FALSE]\n"),

 	make_option(c("--fithicpath"), type="character", default=NA, help="Path to the fithic code [Must for fithic section]
 		[Check for details https://github.com/ay-lab/fithic]\n"),
 	
 	make_option(c("--pythonpath"), type="character", default=NA, help="Path to the python (>= 3 version) [Must for fithic section]\n"),

 	make_option(c("--sthread"), type="integer", default=1, help="Number of threads to be used for parallel sample processing. Can be used using cis/trans/fithic section 
 		[default 1]\n"),
 	
 	make_option(c("--cthread"), type="integer", default=1, help="Number of threads to be used for parallel chromosome processing per sample Can be used using cis/trans section 
 		[default 1]\n"),
 	
 	make_option(c("--pthread"), type="integer", default=1, help="Number of threads to be used for pca calculation per chromosome per sample Can be used using cis/trans section 
 		[default 1]\n
 	Note: The total number of threads will be used is = sthread X cthread X pthread\n")
)
opt <- parse_args(OptionParser(option_list=option_list))

inputfile  <- as.character(opt$file)
chrfile    <- as.character(opt$chrlist)
rowmrge    <- as.integer(opt$rowmrge)
pcatype    <- as.character(opt$pcatype)
genome     <- as.character(opt$genome)
sthread    <- as.integer(opt$sthread)
cthread    <- as.integer(opt$cthread)
pthread    <- as.integer(opt$pthread)
ebackgrnd  <- as.integer(opt$ebackgrnd)
rzscore    <- as.numeric(opt$rzsc)
szscore    <- as.numeric(opt$szsc)
fdr_thr    <- as.numeric(opt$fdr)
gscore     <- as.numeric(opt$gscore)
diffdir    <- as.character(opt$diffdir)
subnum     <- as.integer(opt$subnum)
count_thr  <- as.integer(opt$minc)
dist_thr   <- as.integer(opt$maxd)
distclust  <- as.integer(opt$distclust)
numberclust<- as.integer(opt$numberclust)
bdgfile    <- as.character(opt$bdgfile)
gfolder	   <- as.character(opt$gfolder)
fithicpath <- as.character(opt$fithicpath)
pythonpath <- as.character(opt$pythonpath)
region 	   <- as.character(opt$region)
pcgrp 	   <- as.character(opt$pcgroup)
interaction<- as.character(opt$interaction)
cells	   <- as.character(opt$cells)
dirovwt    <- opt$dirovwt
exclA      <- opt$exclA
pcscore    <- opt$pcscore
compare    <- opt$compare

pc <- as.integer(opt$pc)

if (sthread > 1) {
	cl_sthread<- parallel::makeCluster(sthread)
	parallel::clusterExport(cl_sthread, c("expectedInteraction","ijk2matfunc_cis","extractTrans","ijk2matfunc_trans","mat2fbm"))
}

#Read input file
data <- read.table(paste0(inputfile), h=F, as.is=T)
colnames(data) <- c("mat","bed","prefix","prefix.master")
prefix_master  <- unique(data$prefix.master)

if (pcatype == "cis" | pcatype == "both") {
	if (sthread > 1) {
		invisible(parallel::parLapply(cl_sthread, 1:nrow(data), readfilesintra, data, pc, eigendiv, ebackgrnd, cthread, pthread, dirovwt))
	} else {
		invisible(lapply(1:nrow(data), readfilesintra, data, pc, ebackgrnd, cthread, pthread, dirovwt))
	}
} else if (pcatype == "trans" | pcatype == "both") {
	if (sthread > 1) {
		invisible(parallel::parLapply(cl_sthread, 1:nrow(data), readfilesinter, data, pc, eigendiv, cthread, pthread, rowmrge, dirovwt))
	} else {
		lapply(1:nrow(data), readfilesinter, data, pc, cthread, pthread, rowmrge, dirovwt)
	}
}

#Select and flip the PC coordinates
if (pcatype == "select") {
	if(is.na(genome)) {
		stop("Please provide a valid genome id")
	} else {
		pcselect(data, genome, pc, dirovwt, gfolder)
	}
}

if (pcatype == "analyze") {

	cat ("Running intra sample differential calls for ",prefix_master," samples\n")
	pcanalyze(data, dirovwt, diffdir, rzscore, szscore, fdr_thr, distclust, numberclust)
}

if (pcatype == "subcomp") {

	cat ("Finding sub-compartments from intra sample differential calls for ",prefix_master," samples\n")
	subcompartment(data, diffdir, subnum)
}

if (pcatype == "fithic") {

	cat ("Finding significant loops from intra sample ",prefix_master," replicates\n")
	if (is.na(fithicpath)) {
		stop("Please provide a valid fithic path")
	} else if (is.na(pythonpath)) {
		stop("Please provide a valid python (>= 3 version) path")
	} else {
		if (sthread > 1) {
			fithicformat(data, diffdir, fithicpath, pythonpath, fdr_thr, dist_thr, dirovwt, cl_sthread)
		} else {
			fithicformat(data, diffdir, fithicpath, pythonpath, fdr_thr, dist_thr, dirovwt, cl_sthread=NA)
		}
	}	
}

if (pcatype == "dloop") {
	
	cat ("Finding differential loops across ",prefix_master," samples anchored at least one side within differential compartments\n")
	compartmentLoop(data, diffdir, fdr_thr, count_thr, dist_thr)
}

if (pcatype == "viz") {
	
	cat ("Creating intra sample vizualization for ",prefix_master," samples\n")
	if(is.na(genome)) {
		stop("Please provide a valid genome id")
	} else {
		generateTrackfiles(data, diffdir, genome, bdgfile, pcgrp, fdr_thr, numberclust)
	}
}

if (pcatype == "enrich") {
	
	if(is.na(genome)) {
		stop("Please provide a valid genome id")
	} else {
		cat ("Finding gene enrichment for ",prefix_master," samples\n")
		if (compare) {
			if (cells == "X,Y") {
				stop("Provide a valid sample name\n")
			}
		}
		geneEnrichment(data, diffdir, genome, exclA, region, pcgrp, interaction, pcscore, compare, cells, fdr_thr)
	}
}


if (sthread > 1) {
	parallel::stopCluster(cl_sthread)
}
