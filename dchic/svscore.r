upperqunatmean <- function(x) {
  #mean(x[x >= quantile(x, 0.99)])
  y <- as.vector(quantile(x, c(89:99)/100)) 
  z <- 1:length(y)
  return(lm(y ~ z)$coef[[2]])
}

svscore <- function(mat, bed, prefix, ice) {
 
  folder <- paste0(prefix,"_svscore")
  if (!file.exists(folder)) {
    dir.create(folder)
  }
  setwd(folder)
  ice_command <- paste0(ice," --results_filename ",prefix,".iced --max_iter 250 --filter_low_counts_perc 0.02 --filter_high_counts_perc 0 --output-bias 1 --verbose 1 ",mat)

  if (!file.exists(paste0(prefix,".iced.biases"))) {
    cat ("Running",ice_command,"\n")
    system(ice_command, wait=T)
  }
  
  bed <- read.table(bed, h=F, as.is=T)
  colnames(bed) <- c("chr","start","end","index")
  rownames(bed) <- paste0(bed$chr,"_",bed$start)
  bed[,"bias"]  <- 0
  bed[,"zoe"]   <- 0

  bias <- read.table(paste0(prefix,".iced.biases"), h=F, as.is=T)
  bed$bias[1:nrow(bias)] <- bias$V1
  bed[bed$bias == "NaN", "bias"] <- 0
 
  mat <- as.data.frame(data.table::fread(mat, h=F), as.is=T)
  colnames(mat) <- c("indexA","indexB","count")
  mat[,"chrA"]  <- bed$chr[c(mat$indexA)]
  mat[,"posA"]  <- bed$start[c(mat$indexA)]
  mat[,"biasA"] <- bed$bias[c(mat$indexA)]
  mat[,"chrB"]  <- bed$chr[c(mat$indexB)]
  mat[,"posB"]  <- bed$start[c(mat$indexB)]
  mat[,"biasB"] <- bed$bias[c(mat$indexB)]
  mat <- mat[mat$chrA == mat$chrB,]
  mat[,"dist"] <- abs(mat$posA - mat$posB)
  data <- mat[mat$biasA >= 0.9 & mat$biasA <= 1 & mat$biasB >= 0.9 & mat$biasB <= 1, c("count","dist")] 
  data <- data[order(data$dist),]
  obj.sm <- smooth.spline(x=data$dist, y=data$count)
  mat[,"expect"] <- predict(obj.sm, mat$dist)$y
  mat[,"OE"]  <- mat$count/mat$expect
  oe.vec  <- mat[mat$biasA >= 0.9 & mat$biasA <= 1 & mat$biasB >= 0.9 & mat$biasB <= 1,]$OE
  oe.mean <- mean(oe.vec)
  oe.sd   <- sd(oe.vec)
  mat[,"ZOE"] <- (mat$OE - oe.mean)/oe.sd
  zoe.data  <- data.frame(pos=c(mat$posA,mat$posB), zoe=c(mat$ZOE,mat$ZOE), dist=c(mat$dist,mat$dist), chr=c(mat$chrA,mat$chrB))
  zoe.value <- aggregate(zoe ~ chr + pos, upperqunatmean, data=zoe.data)
  bed[paste0(zoe.value$chr,"_",zoe.value$pos), "zoe"] <- zoe.value$zoe
  setwd("../")
  return(bed)

}  

badregion <- function(data, top, prefix) {

  chrom.name <- unique(data$chr)
  print (chrom.name)
  chrom.list <- list()
  for (i in 1:length(chrom.name)) {
    chrom.list[[i]] <- data[data$chr == as.character(chrom.name[i]),]
    chrom.list[[i]] <- chrom.list[[i]][order(-chrom.list[[i]]$zoe),]
    chrom.list[[i]][,"top"] <- rep(c(1:100),each=ceiling(nrow(chrom.list[[i]])*top))[1:nrow(chrom.list[[i]])]
    #print (chrom.list[[i]])
  }
  chrom.list <- as.data.frame(do.call(rbind, chrom.list), as.is=T)
  chrom.aggr <- aggregate(zoe ~ chr + top, mean, data=chrom.list)
  chrom.aggr$zoe <- round(chrom.aggr$zoe, 3)

  chrom.mat <- matrix(0,100,length(chrom.name))
  rownames(chrom.mat) <- 1:100
  colnames(chrom.mat) <- chrom.name
  for(i in 1:length(chrom.name)) {
    df <- chrom.aggr[chrom.aggr$chr == as.character(chrom.name[i]),]
    df <- df[order(df$top),]
    chrom.mat[1:max(df$top),i] <- df$zoe
  }
  data <- data.frame(top1perc=chrom.mat[1,])
  colnames(data) <- paste0(prefix)
  return(data)
}

runsvscore <- function(file, ice, projname) {
 
  sample.file <- read.table(file, h=T, as.is=T)
  sample.data <- list()
  for(i in 1:nrow(sample.file)) {
    print (sample.file$prefix[i])
    mat <- normalizePath(sample.file$mat[i])
    bed <- normalizePath(sample.file$bed[i])
    prefix <- sample.file$prefix[i]
    data <- svscore(mat, bed, prefix, ice)
    sample.data[[i]] <- badregion(data, top=0.01, prefix)
  }
  sample.data <- do.call(cbind, sample.data)
  outputfile <- paste0(projname, ".svscore.txt")
  write.table(sample.data, file=outputfile, sep="\t", quote=F)

}

## Need ice python script path defined
## and .txt file. 
ice <- normalizePath("ice")

args = commandArgs(trailingOnly=TRUE) 
inputfile = args[length(args)-1]
projname = args[length(args)]

runsvscore(inputfile, ice, projname)

