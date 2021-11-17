## Total interaction
expectedInteraction <- function(i, dist_param, bed_nrow, start_resolution) {
 
  vecInt <- bed_nrow$nrow - (dist_param$dist[i]/start_resolution)
  vecInt <- sum(vecInt[vecInt > 0])
  return(vecInt)
}

## This function converts the Hi-C txt file to matrix file
ijk2matfunc <- function(i, chr, pfx) {

   df  <- data.table::fread(paste0(pfx,"_mat/",chr[i],".txt"), h=T)
   df  <- as.data.frame(df)
   bed <- read.table(paste0(pfx,"_mat/",chr[i],".bed"), h=T)
   mat <- ijk2mat(as.matrix(df[,c("A","B","WeightOE")]), nrow(bed), nrow(bed))
   mat <- as.data.frame(mat)
   colnames(mat) <- paste0(bed$chr,"-",bed$start)
   rownames(mat) <- paste0(bed$chr,"-",bed$start) 
   mat.cor <- suppressWarnings(cor(mat))
   mat.cor <- as.data.frame(mat.cor)
   mat.cor[is.na(mat.cor)] <- 0
   data.table::fwrite(mat.cor, file=paste0(pfx,"_mat/",chr[i],".matrix"), row.names=T, col.names=T, sep="\t", quote=F)
   cat ("Wrote ",pfx,"_mat/",chr[i],".matrix file\n")
}

convertmatrix <- function(prefix) {

  files <- list.files(paste0(prefix,"_mat"), pattern=".txt")
  files <- gsub(".txt", "", files)

  cat ("Found the following chromosomes\n")
  cat (files,"\n")
  invisible(lapply(c(1:length(files)), ijk2matfunc, files, prefix))
  cat ("convert matrix part finished\n")
}

## Read the matrix and bed file
readfiles <- function(mat, bed, prefix) {

  mat <- normalizePath(mat)
  bed <- normalizePath(bed)
  cat ("Reading ",mat," Hi-C matrix file\n")
  cat ("Reading ",bed," Hi-C bed file\n")
  mat <- data.table::fread(mat, h=F)
  mat <- as.data.frame(mat)
  bed <- read.table(bed, h=F)
  colnames(mat) <- c("A","B","Weight")
  colnames(bed) <- c("chr","start","end","index")
  chr.hash <- hashmap::hashmap(bed$index, as.vector(bed$chr))
  pos.hash <- hashmap::hashmap(bed$index, c(bed$start))
  mat[,"chr1"] <- chr.hash[[mat$A]]
  mat[,"pos1"] <- pos.hash[[mat$A]]
  mat[,"chr2"] <- chr.hash[[mat$B]]
  mat[,"pos2"] <- pos.hash[[mat$B]]
  mat <- mat[mat$chr1==mat$chr2,]
  mat[,"dist"] <- abs(mat$pos1 - mat$pos2) 
  print (head(mat))

  dist_param <- aggregate(Weight ~ dist, sum, data=mat)
  print (head(dist_param))  
  
  bed_nrow <- bed
  bed_nrow[,"nrow"] <- 1
  bed_nrow <- aggregate(nrow ~ chr, sum, data=bed_nrow)
  expInt <- lapply(c(1:nrow(dist_param)), expectedInteraction, dist_param, bed_nrow, (bed$end[1] - bed$start[1]))  
  expInt <- do.call(rbind, expInt)
  dist_param[,"totalpairinteraction"] <- expInt
  dist_param[,"expcc"] <- dist_param$Weight/dist_param$totalpairinteraction

  if (!file.exists(paste0(prefix,"_mat"))) {
    dir.create(paste0(prefix,"_mat"))
  }
 
  #write.table(dist_param, file=paste0(prefix,"_mat/ExpectedInteraction.txt"), row.names=F, sep="\t", quote=F)
  expcc_hash <- hashmap::hashmap(dist_param$dist, dist_param$expcc)
  mat[,"WeightOE"] <- mat$Weight/expcc_hash[[mat$dist]]

  chr.uniq <- unique(bed$chr)
  for(i in 1:length(chr.uniq)) {
    cat ("Writing ",as.character(chr.uniq[i]),".txt file\n")
    b <- bed[bed$chr==chr.uniq[i],]
    b[,"index"] <- 1:nrow(b)
    index.hash  <- hashmap::hashmap(b$start, b$index)
    m <- mat[(mat$chr1==chr.uniq[i] & mat$chr2==chr.uniq[i]), c(4,5,6,7,9,1,2)]
    m[,"A"] <- index.hash[[m$pos1]]
    m[,"B"] <- index.hash[[m$pos2]] 
    data.table::fwrite(m, file=paste0(prefix,"_mat/",chr.uniq[i],".txt"), row.names=F, col.names=T, sep="\t", quote=F)
    data.table::fwrite(b, file=paste0(prefix,"_mat/",chr.uniq[i],".bed"), row.names=F, col.names=T, sep="\t", quote=F)
  }
}


#mat <- c("../Mouse_Result/Matrix/mESC_Bonev_40Kb.matrix","../Mouse_Result/Matrix/NPC_Bonev_40Kb.matrix","../Mouse_Result/Matrix/CN_Bonev_40Kb.matrix")
#bed <- c("../Mouse_Result/Matrix/mESC_Bonev_40Kb_abs.bed","../Mouse_Result/Matrix/NPC_Bonev_40Kb_abs.bed","../Mouse_Result/Matrix/CN_Bonev_40Kb_abs.bed")
#prefix <- c("mESC","NPC","CN")

hicmat <- function(mat, bed, prefix) {
  for(n in 1:length(prefix)){
    cat ("Running ",prefix[n],"\n") 
    readfiles(mat[n], bed[n], prefix[n])
    convertmatrix(prefix[n])
  }
}
