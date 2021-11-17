## This function helps to find any executable
.findExecutable <- function(exe, interactive=TRUE) {
  path <- Sys.which(exe)
  if(all(path=="")) {
    if(interactive) stop("Executable for ", paste(exe, collapse=" or "), " not found! Please make sure that the software is correctly installed and, if necessary, path variables are set.", call.=FALSE)
    return(character(0))
  }
  path[which(path!="")[1]]
}

## Cluster Function - Find compartments enriched in a certain experiment
clust <- function(folder, chr, orient="Active2Inactive") {

  bdg <- list.files(paste0(folder,"/",chr), pattern="HMFA_chr_")
  pc  <- list()
  for(i in 1:length(bdg)) {
    if (i == 1) { 
      bed <- read.table(paste0(folder,"/",chr,"/",bdg[i]), h=F, as.is=T)
      pc[[i]] <- round(bed[,4],3)
      bed <- bed[,-4] 
      colnames(bed) <- c("chr","start","end")
    } else {
      pc[[i]] <- round(read.table(paste0(folder,"/",chr,"/",bdg[i]), h=F, as.is=T)[,4],3)
    }
  }
  pc  <- as.data.frame(do.call(cbind, pc))
  colnames(pc) <- do.call(rbind, strsplit(bdg,"_|[.]"))[,5]
  magnitude <- list()
  namelist  <- list()
  compchng  <- list()
  signchg   <- list()
  pcdata    <- list()
  for(i in 1:nrow(pc)) {
    v <- data.frame(value=as.numeric(pc[i,]),name=colnames(pc[i,]))
    if (orient == "Inactive2Active") {
      v <- v[order(v$value),]
      s <- sign(v$value)
      s[s==1]  <- "P"
      s[s==-1] <- "N"
      magnitude[[i]] <- round((mean(dist(as.numeric(v$value)))),3)
      namelist[[i]]  <- paste0(v$name, collapse="->")
      signchg[[i]]   <- paste0(s,collapse="->")
      pcdata[[i]]    <- as.numeric(v$value)
      if (v$value[length(v$value)] > 0 & v$value[1] > 0) {
        compchng[[i]] <- "L.A->H.A"
      } else if (v$value[length(v$value)] <= 0 & v$value[1] < 0) {
        compchng[[i]] <- "L.B->H.B"
      } else if (v$value[length(v$value)] > 0 & v$value[1] <= 0) {
        compchng[[i]] <- "B->A"
      } else {
        compchng[[i]] <- "Undefined"
      }
    } else if (orient == "Active2Inactive") {
      v <- v[order(-v$value),]
      s <- sign(v$value)
      s[s==1]  <- "P"
      s[s==-1] <- "N"
      magnitude[[i]] <- round((mean(dist(as.numeric(v$value)))),3)
      namelist[[i]]  <- paste0(v$name, collapse="->")
      signchg[[i]]   <- paste0(s,collapse="->")
      pcdata[[i]]    <- as.numeric(v$value)
      if (v$value[length(v$value)] > 0 & v$value[1] > 0) {
        compchng[[i]] <- "H.A->L.A"
      } else if (v$value[length(v$value)] < 0 & v$value[1] <= 0) {
        compchng[[i]] <- "H.B->L.B"
      } else if (v$value[length(v$value)] <= 0 & v$value[1] > 0) {
        compchng[[i]] <- "A->B"
      } else {
        compchng[[i]] <- "Undefined"
      }
    }
  }
  magnitude <- unlist(magnitude)
  namelist  <- unlist(namelist)
  compchng  <- unlist(compchng)
  signchg   <- unlist(signchg)
  pcdata    <- as.data.frame(do.call(rbind, pcdata))
  colnames(pcdata) <- paste0("S",c(1:ncol(pcdata)))
  pcdata[,"magnitude"] <- magnitude
  pcdata[,"namelist"]  <- namelist
  pcdata[,"compchng"]  <- compchng
  pcdata[,"signlist"]  <- signchg
  pcdata <- data.frame(bed, pcdata)
  write.table(pcdata, file=paste0(chr,".pc_change.txt"), row.names=F, sep="\t", quote=F)
}

separation <- function(chr, include, exclude, type="both", orient="Active2Inactive") {

  pcdata   <- read.table(paste0(chr,".pc_change.txt"), h=T, as.is=T)
  namelist <- pcdata$namelist
  signlist <- pcdata$signlist
  namelist <- do.call(rbind,strsplit(namelist,"->"))
  signlist <- do.call(rbind,strsplit(signlist,"->"))
  criticalcompt <- list()
  k <- 1
  for(i in 1:nrow(namelist)) {
    w <- length(namelist[i,]):1
    s <- signlist[i,]
    s[s=="P"] <- 1 
    s[s=="N"] <- -1
    names(w) <- namelist[i,]
    names(s) <- namelist[i,]
    if (type == "direction" & (min(w[include]) > max(w[exclude]))) {
      criticalcompt[[k]] <- pcdata[i,]
      k <- k + 1
    } else if (type == "sign" & min(s[include]) == 1 & orient == "Active2Inactive") {
      criticalcompt[[k]] <- pcdata[i,]
      k <- k + 1
    } else if (type == "sign" & min(s[include]) == -1 & orient == "Inactive2Active") {
      criticalcompt[[k]] <- pcdata[i,]
      k <- k + 1
    } else if (type == "both" & ((min(w[include]) > max(w[exclude])) & min(s[include]) == 1 & orient == "Active2Inactive")) {
      criticalcompt[[k]] <- pcdata[i,]
      k <- k + 1
    } else if (type == "both" & ((min(w[include]) > max(w[exclude])) & min(s[include]) == -1 & orient == "Inactive2Active")) {
      criticalcompt[[k]] <- pcdata[i,]
      k <- k + 1
    }
  }
  criticalcompt <- do.call(rbind, criticalcompt)
  return(criticalcompt)
}

mapGenes <- function(dcomp.file, gene_pos, slack) {
  tmp <- tempfile(pattern = "temp", tmpdir = tempdir(), fileext = ".bed")
  dcomp<- read.table(dcomp.file, h=T, comment.char="", as.is=T)
  #dcomp<- dcomp[dcomp$padj < padj,] # & dcomp$dZsc > dzsc
  colnames(dcomp)[1] <- "chr"
  dcomp <- dcomp[order(dcomp$chr),]
  slackwin <- (dcomp$end[1] - dcomp$start[1]) * slack
  slackwin <- as.integer(slackwin)
  rownames(dcomp) <- paste0(dcomp$chr,":",dcomp$start,"-",dcomp$end)
  write.table(dcomp[,1:3], file=tmp, col.names=F, row.names=F, sep="\t", quote=F)
  #cmd <- paste0(.findExecutable("bedtools")," map -a ",tmp," -b ",gene_pos," -c 4 -o collapse -null NA")
  cmd <- paste0(.findExecutable("bedtools")," window -a ",tmp," -b ",gene_pos," -w ",slackwin," |",.findExecutable("bedtools")," groupby -g 1,2,3 -c 7 -o distinct")
  df <- read.table(text=system(cmd, wait=T,intern=T),h=F)
  rownames(df) <- paste0(df$V1,":",df$V2,"-",df$V3)
  dcomp[,"gene"] <- NA
  dcomp[,"gene"] <- df[rownames(dcomp),4]
  write.table(na.omit(dcomp), file=paste0("Genes.",basename(dcomp.file)), row.names=F, sep="\t", quote=F)
}

args = commandArgs(trailingOnly=TRUE) 
cluster <- args[length(args)]
print(cluster)

if (cluster == 1) {
  folder <- args[length(args)-1]
  include <- strsplit(args[length(args)-2], ",")[[1]]
  exclude <- strsplit(args[length(args)-3], ",")[[1]]
  orientation <- args[length(args)-4]
  typeOption <- args[length(args)-5]
  chrlist <- list.files(folder, pattern="chr_")
  chrlist <- chrlist[which(chrlist != "chr_info.txt")]
  crtcomp <- list()
  for(j in 1:length(chrlist)) {
    data <- clust(folder, chrlist[j], orient=orientation)
    crtcomp[[j]] <- separation(chrlist[j], include, exclude, type=typeOption, orient=orientation)
  }
  crtcomp <- do.call(rbind, crtcomp)
  #print (crtcomp)
  write.table(crtcomp, file="tmp", row.names=F, sep="\t", quote=F)
} else {
  dcompFile <- args[length(args)-1]
  genePosFile <- args[length(args)-2]
  slack <- as.integer(args[length(args)-3])
  mapGenes(dcompFile, genePosFile, slack)
}
  
