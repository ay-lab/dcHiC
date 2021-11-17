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

## Process the downloaded data
process <- function(genome, resolution, folder=NA) {

  if (is.na(folder)) {
    folder <- paste0(genome,"_goldenpathData")
    if (!file.exists(folder)) {
      dir.create(folder)
    }
    folder <- normalizePath(folder)
  } else {
    folder <- normalizePath(folder)
  }

  current_path <- getwd()
  setwd(folder)
  datadownload(genome)
  setwd(current_path)

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
}

## Perform chromosome wise regression to determine best pc and sign
pcselect <- function(obj, chr, genome, resolution, folder=NA, lnkeep=0.5, tskeep=0.1, gckeep=0.1, component="ind") {

  if (is.na(folder)) {
    folder <- paste0(genome,"_goldenpathData_",as.numeric(resolution))
    if (!file.exists(folder)) {
      dir.create(folder)
    }
    folder <- normalizePath(folder)
  } else {
    folder <- normalizePath(folder)
  }

  ## Process the files
  process(genome, resolution, folder)

  gcpt_master <- read.table(paste0(folder,"/",genome,".GCpt.tss.bedGraph"), h=F, as.is=T)
  colnames(gcpt_master) <- c("chr","start","end","gc","tss")
  gcpt_master <- gcpt_master[gcpt_master$chr == chr,]
  rownames(gcpt_master) <- paste0(gcpt_master$chr,"-",gcpt_master$start)

  ## When component is partial
  if (component == "partial") {
    ## Number of samples
    nsample <- ncol(obj$partial[[1]][,1,])
    npc <- ncol(obj$partial[[1]][,,1])
    weight <- as.vector(round(obj$eig[1:npc,2],3)/sum(round(obj$eig[1:npc,2],3)))
    gc  <- list()
    ln  <- list()
    ts  <- list()
    for (i in 1:nsample) {
      pc <- list()
      for (j in 1:npc) {
        pc[[j]] <- as.vector(obj$partial[[1]][,j,i])
      }
      pc <- as.data.frame(do.call(cbind, pc))
      colnames(pc) <- paste0("pc",c(1:npc))
      gcpc<- data.frame(gcpt_master[names(obj$partial[[1]][,1,i]),], pc)
      ln[[i]] <- data.frame(t(cor(cbind(len=1:nrow(gcpc),gcpc[,6:ncol(gcpc)]),method="spearman")[,1]),id=i,chr=chr)
      gc[[i]]  <- data.frame(t(cor(gcpc[,c(4,6:ncol(gcpc))], method = 'spearman')[,1]),id=i,chr=chr)
      ts[[i]] <- data.frame(t(cor(gcpc[,c(5:ncol(gcpc))], method = 'spearman')[,1]),id=i,chr=chr)
    }
    ln <- as.data.frame(do.call(rbind, ln))
    gc <- as.data.frame(do.call(rbind, gc))
    ts <- as.data.frame(do.call(rbind, ts))
    ln[,"sample"] <- rownames(obj$group$coord[[1]])
    gc[,"sample"] <- rownames(obj$group$coord[[1]])
    ts[,"sample"] <- rownames(obj$group$coord[[1]])
 
    #print (gc)
    #print (ts) 
    pc.select <- list()
    for (i in 1:nrow(ln)) { 
      m <- data.frame(t(abs(ln[i,2:(ncol(ln)-3)])), t(abs(gc[i,2:(ncol(gc)-3)])), t(ts[i,2:(ncol(ts)-3)]))
      rownames(m)  <- 1:nrow(m)
      colnames(m)  <- c("ln","gc","tss")
      m[,"gcsign"] <- t(sign(gc[i,2:(ncol(gc)-3)]))[,1]
      m$tss <- m$tss * m$gcsign
      m[,"pc"]     <- 1:nrow(m)
      m[,"score"]  <- (((1 - m$ln) + m$gc + m$tss)/3) + weight
      m[,"lnkeep"] <- "yes"
      m[,"tskeep"] <- "no"
      m[,"gckeep"] <- "yes"
      m[,"sample"] <- rownames(obj$group$coord[[1]])[i]
      while (nrow(m[m$ln <= lnkeep,]) == 0) {
        cat ("All the PC's are correlated above the lnkeep threshold ",lnkeep,". Increasing the lnkeep threshold. PC's are highly correlated with chromsome arm/length.\n")
        lnkeep <- lnkeep + 0.05
      }
      m[m$ln > lnkeep, "lnkeep"] <- "no"
    
      while (nrow(m[m$tss > tskeep,]) == 0) {
        cat ("All the PC's are weekly correlated with transcription start site. PC's should have high correlation with tss. Lowering tskeep.\n")
        tskeep <- tskeep - 0.05
      }
      m[m$tss > tskeep, "tskeep"] <- "yes"    

      while (nrow(m[m$gc > gckeep,]) == 0) {
        cat ("All the PC's are weekly correlated with gc content. PC's should have high correlation with gc content. Lowering gckeep.\n")
        gckeep <- gckeep - 0.05
      }
      m[m$gc > gckeep, "gckeep"] <- "yes"
      m[,"pcselect"] <- "no"

      while (nrow(m[m$lnkeep == "yes" & m$tskeep == "yes" & m$gckeep == "yes",]) == 0) {
        cat ("Increasing the lnkeep threshold.\n")
        lnkeep <- lnkeep + 0.05 
        m[m$ln < lnkeep, "lnkeep"] <- "yes"
      }

      print (i)
      print (m)
      cat ("\n")

      pc.select[[i]] <- m 
      m <- m[m$lnkeep == "yes" & (m$tskeep == "yes" | m$gckeep == "yes"),]
      m <- m[order(-m$score),]
      pc.select[[i]][pc.select[[i]]$pc==m$pc[1],"pcselect"] <- "yes"
    }
    pc.select <- as.data.frame(do.call(rbind, pc.select))
    write.table(pc.select, file="PCselection.txt", row.names=F, col.names=T, sep="\t", quote=F)
    pc.select <- pc.select[pc.select$pcselect == "yes",]  

    pc <- list()
    for(i in 1:nrow(pc.select)) {
      pc[[i]] <- data.frame(gcpt_master[names(obj$partial[[1]][,1,i]),],pc=(obj$partial[[1]][,pc.select$pc[i],i] * pc.select$gcsign[i]))
      #pc[[i]] <- (obj$partial[[1]][,pc.select$pc[i],i] * pc.select$gcsign[i])
    }
    #print (cor(do.call(cbind,pc)))
    #return(as.data.frame(do.call(cbind,pc), as.is=T))
    #return(pc)
  } else if (component == "ind") {
    npc <- ncol(obj$ind[[1]])
    weight <- as.vector(round(obj$eig[1:npc,2],3)/sum(round(obj$eig[1:npc,2],3)))
    gc  <- list()
    ln  <- list()
    ts  <- list() 
    pc  <- list()
    for (i in 1:npc) {
      pc[[i]] <- as.vector(obj$ind[[1]][,i])
    }
    pc <- as.data.frame(do.call(cbind, pc))
    colnames(pc) <- paste0("pc",c(1:npc))
    gcpc<- data.frame(gcpt_master[names(obj$ind[[1]][,1]),], pc)
    ln[[1]] <- data.frame(t(cor(cbind(len=1:nrow(gcpc),gcpc[,6:ncol(gcpc)]),method="spearman")[,1]),id=i,chr=chr)
    gc[[1]]  <- data.frame(t(cor(gcpc[,c(4,6:ncol(gcpc))], method = 'spearman')[,1]),id=i,chr=chr)
    ts[[1]] <- data.frame(t(cor(gcpc[,c(5:ncol(gcpc))], method = 'spearman')[,1]),id=i,chr=chr)
    ln <- as.data.frame(do.call(rbind, ln))
    gc <- as.data.frame(do.call(rbind, gc))
    ts <- as.data.frame(do.call(rbind, ts))
    ln[,"sample"] <- "Average"
    gc[,"sample"] <- "Average"
    ts[,"sample"] <- "Average"
    pc.select <- list()
    for (i in 1:nrow(ln)) {
      m <- data.frame(t(abs(ln[i,2:(ncol(ln)-3)])), t(abs(gc[i,2:(ncol(gc)-3)])), t(ts[i,2:(ncol(ts)-3)]))
      rownames(m)  <- 1:nrow(m)
      colnames(m)  <- c("ln","gc","tss")
      m[,"gcsign"] <- t(sign(gc[i,2:(ncol(gc)-3)]))[,1]
      m$tss <- m$tss * m$gcsign
      m[,"pc"]     <- 1:nrow(m)
      m[,"score"]  <- (((1 - m$ln) + m$gc + m$tss)/3) + weight
      m[,"lnkeep"] <- "yes"
      m[,"tskeep"] <- "no"
      m[,"gckeep"] <- "yes"
      m[,"sample"] <- "Average"
      while (nrow(m[m$ln <= lnkeep,]) == 0) {
        cat ("All the PC's are correlated above the lnkeep threshold ",lnkeep,". Increasing the lnkeep threshold. PC's are highly correlated with chromsome arm/length.\n")
        lnkeep <- lnkeep + 0.05
      }
      m[m$ln > lnkeep, "lnkeep"] <- "no"

      while (nrow(m[m$tss > tskeep,]) == 0) {
        cat ("All the PC's are weekly correlated with transcription start site. PC's should have high correlation with tss. Lowering tskeep.\n")
        tskeep <- tskeep - 0.05
      }
      m[m$tss > tskeep, "tskeep"] <- "yes"

      while (nrow(m[m$gc > gckeep,]) == 0) {
        cat ("All the PC's are weekly correlated with gc content. PC's should have high correlation with gc content. Lowering gckeep.\n")
        gckeep <- gckeep - 0.05
      }
      m[m$gc > gckeep, "gckeep"] <- "yes"
      m[,"pcselect"] <- "no"

      while (nrow(m[m$lnkeep == "yes" & m$tskeep == "yes" & m$gckeep == "yes",]) == 0) {
        cat ("Increasing the lnkeep threshold.\n")
        lnkeep <- lnkeep + 0.05
        m[m$ln < lnkeep, "lnkeep"] <- "yes"
      }

      print (i)
      print (m)
      cat ("\n")

      pc.select[[i]] <- m
      m <- m[m$lnkeep == "yes" & (m$tskeep == "yes" | m$gckeep == "yes"),]
      m <- m[order(-m$score),]
      pc.select[[i]][pc.select[[i]]$pc==m$pc[1],"pcselect"] <- "yes"
    }

    pc.select <- as.data.frame(do.call(rbind, pc.select))
    write.table(pc.select, file="PCselection.txt", row.names=F, col.names=T, sep="\t", quote=F)
    pc.select <- pc.select[pc.select$pcselect == "yes",]

    nsample <- ncol(obj$partial[[1]][,1,])
    pc <- list()
    for(i in 1:nsample) {
      pc[[i]] <- data.frame(gcpt_master[names(obj$partial[[1]][,1,i]),],pc=(obj$partial[[1]][,pc.select$pc[1],i] * pc.select$gcsign[1]))
      #pc[[i]] <- (obj$partial[[1]][,pc.select$pc[1],i] * pc.select$gcsign[1])
    }
    #pc <- as.data.frame(do.call(cbind,pc), as.is=T)
    #colnames(pc) <- rownames(obj$group$coord[[1]])
    #write.table(pc, file="PCs.txt", sep="\t", quote=F)
    #return(as.data.frame(do.call(cbind,pc), as.is=T))
    return(pc)
  }
} 

## obj : HMFA object file
## genome : hg38/hg19/mm9/mm10
## chr : chromosome name for which HMFA was ran
## resolution : HiC resolution
## folder : path to the already existing _goldenpath folder generated by armcorrection. If you don't know, keep it as NA
## component : There are two ways to select. First one is to use the partial factor scores "partial" and second one is the average of all partial scores calles "ind" (default")
## This is will return a modfied list object and each item in the list is the sample and can be called as 
## obj <- pcselect(hmfa.obj, "chr22", "hg19", 25000, folder=NA) and the pc values can be accessed like 
## obj[[i]]$pc, where i is the sample number

#obj <- readRDS("/mnt/BioAdHoc/Groups/vd-ay/jeffreywang/DiffComptTesting/EBV_100kb_filtered/SVS_Cutoff_10/chr_3/Rsession_chr3.rds")
#pcselect(obj, "chr3", "hg19", 100000, folder=NA, lnkeep=0.7, component="ind")
#write.table(pcselect(obj, "chr1", "hg19", 1000000, folder=NA, lnkeep=0.6), file="pc.txt", row.names=F, sep="\t", quote=F)
