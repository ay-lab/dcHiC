args = commandArgs(trailingOnly=TRUE) # pcfilenames
path = normalizePath(args[length(args)])
print(path)
source(file.path(path,"createfolder.r"))

distancePointLine <- function(p1x, p1y, slope, intercept) {
  val = abs((p1x * slope * -1) + (p1y) - intercept) / sqrt(slope ** 2 + 1)
  return(val)
}

mdweight <- function(m) {
  d <- list()
  for(i in 1:length(m)) {
    d[[i]] <- sum(abs(m[i] - m))
  }
  d <- unlist(d)
  w <- list()
  for (i in 1:length(d)) {
    w[[i]] <- (length(d) * d[i])/sum(d)
  }
  w <- unlist(w)
  return(w)
}

fitModelChromwise <- function(X,Y,M,S,mcomp=FALSE) {
  
  D <- data.frame(X=(X-mean(X))/sd(X),Y=(Y-mean(Y))/sd(Y))
  model <- lm(Y ~ X, data=D)
  cfs <- coef(model)
  D[,"distance"] <- 0
  for(i in 1:nrow(D)) {
    D$distance[i] <- round(distancePointLine(p1x=D$X[i], p1y=D$Y[i], slope=cfs[2], intercept=cfs[1]),3)
  }
  if (mcomp) {
    return(data.frame(z=c((D$distance - M)/S),d=D$distance))
  } else {
    return(c((D$distance - M)/S))
  }
}

## Input 
## @pcfiles is a vector of chromosome wise pc component file names
## @samplefile is a dataframe with the following format
## mat	bed	replicate	prefix	group
## A_rep1.matrix	A_rep1_abs.bed	A_rep1	A	<optional>
## A_rep2.matrix	A_rep2_abs.bed	A_rep2	A	<optional>
## B_rep1.matrix        B_rep1_abs.bed  B_rep1  B       <optional>
## B_rep2.matrix        B_rep2_abs.bed  B_rep2  B       <optional>
## ..
##
## @resolution is a integer value e.g. 1e5
## @lprmfile is a dataframe with the following format
## m	s	chr
## 0.1	0.2	chr1
## ...
##
## padj is the adjusted pvalue threshold 
## dzsc is the zscore cutoff (deprecated feature)
## mcomp is the multi comparison. Set to FALSE if not required 
#

diffcmp <- function(chr, pcfiles, samplefile, lprmfile, resolution, mcomp=TRUE, padj=1e-2, dzsc=3) {

  createFolder("DifferentialCompartment",1)
  if (length(chr) == 1) {
    createFolder(paste0(chr,"_diffcomp"),0)
  }
  if (file.exists(paste0(lprmfile))) {
    cat ("Learned parameter file found. Using the IHW to boost statistical power\n")
    param <- read.table(lprmfile,h=T, as.is=T)
    combined_data <- list()
    combined_distdf <- list()
    if (ncol(samplefile) == 3) {
      prefix <- unique(samplefile$group)
    } else { 
      prefix <- unique(samplefile$prefix)
    }
    for(i in 1:length(pcfiles)) {
      data <- list()
      l <- 1
      df <- read.table(paste0(pcfiles[i]), h=T, as.is=T)
      for(j in 1:length(prefix)) {
        if (ncol(samplefile) == 3) {
          cat ("Running ",as.character(prefix[j])," ",as.character(pcfiles[i]),"\n")
          selected <- as.vector(samplefile[samplefile$group==prefix[j],]$replicate)
          data[[l]] <- apply(df[,selected], 1, mean)
          #data[[l]] <- apply(df[,c(samplefile[samplefile$group==prefix[j],]$replicate)], 1, mean)
        } else {
          cat ("Running ",as.character(prefix[j])," ",as.character(pcfiles[i]),"\n")
          #print(print(c(samplefile[samplefile$group==prefix[j],]$replicate)))
          if (nrow(samplefile[samplefile$prefix==prefix[j],]) > 1) {
            #print("Entered here.")
            #print(prefix[j])
            #print(samplefile$prefix)
            #print(samplefile[samplefile$prefix==prefix[j],])
            selected <- as.vector(samplefile[samplefile$prefix==prefix[j],]$replicate)
            #print(selected)
            data[[l]] <- apply(df[,selected], 1, mean)
            #head(df[,c(samplefile[samplefile$prefix==prefix[j],]$replicate)])
            #print(c(samplefile[samplefile$prefix==prefix[j],]$replicate))
            #data[[l]] <- apply(df[,c(samplefile[samplefile$prefix==prefix[j],]$replicate)], 1, mean)
          } else { 
            #print("Entered other here.")
            selected <- as.vector(samplefile[samplefile$prefix==prefix[j],]$replicate)
            data[[l]] <- df[,selected]
            #data[[l]] <- df[,c(samplefile[samplefile$prefix==prefix[j],]$replicate)]
          }
        }
        l <- l + 1
      }
      data <- do.call(cbind, data) 
      data <- as.data.frame(data)
      colnames(data) <- prefix
      rownames(data) <- rownames(df) 
      print (head(data))
      distdf <- data
      r1 <- 1
      while (r1 < length(prefix)) {
        r2 <- r1 + 1
        while (r2 <= length(prefix)) {
          cl_md   <- paste0(colnames(data)[r1],"_vs_",colnames(data)[r2],"_md")
          cl_dzsc <- paste0(colnames(data)[r1],"_vs_",colnames(data)[r2],"_dZsc")
          cl_pval <- paste0(colnames(data)[r1],"_vs_",colnames(data)[r2],"_Pval") 
          data[,cl_md]   <- mahalanobis(data[,c(r1,r2)], colMeans(data[,c(r1,r2)]), cov(data[,c(r1,r2)]))
          #data[,cl_md]   <- mahalanobis(scale(data[,c(r1,r2)]), MASS::cov.mcd(scale(data[,c(r1,r2)]))$center, MASS::cov.mcd(scale(data[,c(r1,r2)]))$cov)
          data[,cl_md]   <- data[,cl_md] * mdweight(data[,cl_md])
          data[,cl_pval] <- pchisq(data[,cl_md], 1, lower.tail=F)
          if (mcomp) {
            cl_dist <- paste0(colnames(data)[r1],"_vs_",colnames(data)[r2],"_dist")
            val <- fitModelChromwise(as.vector(data[,r1]),as.vector(data[,r2]),param[param$chr==chr[i],]$m,param[param$chr==chr[i],]$s, mcomp)
            data[,cl_dzsc]   <- val$z
            distdf[,cl_dist] <- val$d
          } else {
            data[,cl_dzsc] <- fitModelChromwise(as.vector(data[,r1]),as.vector(data[,r2]),param[param$chr==chr[i],]$m,param[param$chr==chr[i],]$s, mcomp)
          }
          r2 <- r2 + 1
        }
        r1 <- r1 + 1
      }
      combined_data[[i]] <- data
      #print (ncol(combined_data[[i]]))
      if (mcomp) {
        combined_distdf[[i]] <- distdf
      }
    }
    combined_data <- do.call(rbind, combined_data)
    bed <- as.data.frame(do.call(rbind,strsplit(rownames(combined_data),"-")),stringsAsFactors=F)
    bed$V1 <- as.character(bed$V1)
    bed$V2 <- as.integer(bed$V2)
    colnames(bed) <- c("chr","start")
    bed[,"end"] <- bed$start + as.integer(resolution)
    r1 <- 1
    while (r1 < length(prefix)) {
      r2 <- r1 + 1
      while (r2 <= length(prefix)) {
        cl_md   <- paste0(colnames(combined_data)[r1],"_vs_",colnames(combined_data)[r2],"_md")
        cl_dzsc <- paste0(colnames(combined_data)[r1],"_vs_",colnames(combined_data)[r2],"_dZsc")
        cl_pval <- paste0(colnames(combined_data)[r1],"_vs_",colnames(combined_data)[r2],"_Pval")
        cl_padj <- paste0(colnames(combined_data)[r1],"_vs_",colnames(combined_data)[r2],"_Padj")
        obj <- data.frame(dzsc=combined_data[,cl_dzsc],pval=combined_data[,cl_pval])
        #combined_data[,cl_padj] <- IHW::adj_pvalues(IHW::ihw(pval ~ dzsc, data=obj, alpha=0.01))
        combined_data[,cl_padj] <- invisible(IHW::adj_pvalues(IHW::ihw(pval ~ dzsc, data=obj, alpha=0.01)))
        bed_result <- data.frame(bed,combined_data[,c(colnames(combined_data)[r1],colnames(combined_data)[r2],cl_md,cl_dzsc,cl_pval,cl_padj)])
        bed_result <- bed_result[order(bed_result$chr),]
        colnames(bed_result)[1] <- "#chr"
        colnames(bed_result)[6] <- "mdist"
        colnames(bed_result)[7] <- "dZsc"
        colnames(bed_result)[8] <- "pval"
        colnames(bed_result)[9] <- "padj"
        if (length(chr) == 1) {
          write.table(bed_result, file=paste0(chr,"_diffcomp/",colnames(combined_data)[r1],"_vs_",colnames(combined_data)[r2],"_full_compartment_details.bedGraph"),sep="\t",row.names=F,quote=F)
        } else {
          write.table(bed_result, file=paste0(colnames(combined_data)[r1],"_vs_",colnames(combined_data)[r2],"_full_compartment_details.bedGraph"),sep="\t",row.names=F,quote=F)
        }
        bed_result <- bed_result[bed_result$padj < padj,] # removed dZsc threshold
        if (length(chr) == 1) {
          write.table(bed_result, file=paste0(chr,"_diffcomp/",colnames(combined_data)[r1],"_vs_",colnames(combined_data)[r2],"_differential_compartments.bedGraph"),sep="\t",row.names=F,quote=F)
        } else {
          write.table(bed_result, file=paste0(colnames(combined_data)[r1],"_vs_",colnames(combined_data)[r2],"_differential_compartments.bedGraph"),sep="\t",row.names=F,quote=F)
        }
        #write.table(bed_result, file=paste0(colnames(combined_data)[r1],"_vs_",colnames(combined_data)[r2],"_full_compartment_details.bedGraph"),sep="\t",row.names=F,quote=F)
        #write.table(bed_result, file=paste0(colnames(combined_data)[r1],"_vs_",colnames(combined_data)[r2],"_differential_compartments.bedGraph"),sep="\t",row.names=F,quote=F)
        r2 <- r2 + 1
      }
      r1 <- r1 + 1
    }
    if (mcomp) {
      print("Running MultiComparison")
      combined_distdf <- do.call(rbind, combined_distdf)
      bed <- as.data.frame(do.call(rbind,strsplit(rownames(combined_distdf),"-")),stringsAsFactors=F)
      bed$V1 <- as.character(bed$V1)
      bed$V2 <- as.integer(bed$V2)
      colnames(bed) <- c("chr","start")
      bed[,"end"] <- bed$start + as.integer(resolution)
      #print (head(combined_distdf[,c((1+length(prefix)):ncol(combined_distdf))]))
      distdf <- combined_distdf[,c((1+length(prefix)):ncol(combined_distdf))]
      combined_distdf <- combined_distdf[,-c((1+length(prefix)):ncol(combined_distdf))]
      combined_distdf[,"dist"] <- apply(distdf, 1, mean)
      combined_distdf <- data.frame(bed, combined_distdf)
      combined_distdf[,"mdist"]  <- 0
      #print (head(combined_distdf))
      combined_data <- list()
      for(i in 1:length(chr)) {
        adjstr = paste("chr", chr[i], sep = "")
        data <- combined_distdf[combined_distdf$chr==adjstr,]
        data[,"mdist"]<- mahalanobis(data[,c(4:(ncol(data)-2))], colMeans(data[,c(4:(ncol(data)-2))]), cov(data[,c(4:(ncol(data)-2))]))
        #data[,"mdist"]<- mahalanobis(scale(data[,c(4:(ncol(data)-2))]), MASS::cov.mcd(scale(data[,c(4:(ncol(data)-2))]))$center, MASS::cov.mcd(scale(data[,c(4:(ncol(data)-2))]))$cov)
        data[,"mdist"]<- data[,"mdist"] * mdweight(data[,"mdist"])
        data[,"dZsc"] <- (data$dist - param[param$chr==chr[i],]$m)/param[param$chr==chr[i],]$s
        data[,"pval"] <- pchisq(data[,"mdist"], c(ncol(data)-6), lower.tail=F)    
        combined_data[[i]]  <- data
      }
      combined_data <- do.call(rbind, combined_data)
      obj <- data.frame(dzsc=combined_data[,"dZsc"],pval=combined_data[,"pval"])
      combined_data[,"padj"] <- IHW::adj_pvalues(IHW::ihw(pval ~ dzsc, data=obj, alpha=0.01))
      combined_data <- combined_data[,-c(which(colnames(combined_data)=="dist"))]
      bed_result <- combined_data
      bed_result <- bed_result[order(bed_result$chr),]
      colnames(bed_result)[1] <- "#chr"
      if (length(chr) == 1) {
        write.table(bed_result, file=paste0(chr,"_diffcomp/","MultiComparison_full_compartment_details.bedGraph"),sep="\t",row.names=F,quote=F)
      } else {
        write.table(bed_result, file=paste0("MultiComparison_full_compartment_details.bedGraph"),sep="\t",row.names=F,quote=F)
      }
      bed_result <- bed_result[bed_result$padj < padj,] #removed dZsc threshold
      if (length(chr) == 1) {
        write.table(bed_result, file=paste0(chr,"_diffcomp/","MultiComparison_differential_compartments.bedGraph"),sep="\t",row.names=F,quote=F) 
      } else { 
        write.table(bed_result, file=paste0("MultiComparison_differential_compartments.bedGraph"),sep="\t",row.names=F,quote=F)
      }
      #write.table(bed_result, file=paste0("MultiComparison_full_compartment_details.bedGraph"),sep="\t",row.names=F,quote=F)
      #bed_result <- bed_result[bed_result$padj < padj & bed_result$dZsc > dzsc,]
      #write.table(bed_result, file=paste0("MultiComparison_differential_compartments.bedGraph"),sep="\t",row.names=F,quote=F) 
    }
  } else {
    stop("Please provide a replicate parameter file\n")
  }
}

# Main Method
samplefile = read.table(args[length(args)-1], h=T)
lprmfile = normalizePath(args[length(args)-2])

if (file.exists(paste0(lprmfile))) {
  print("Parameter file found.")
}
res = as.numeric(args[length(args)-3])

chrfile = args[length(args)-4]
if (file.exists(chrfile)) {
  chrvector <- as.character(scan(chrfile, what = character()))
} else {
  chrvector <- c(chrfile) # single chromosome case
}

multicomp = args[length(args)-5]
if (multicomp == 1) {
  mcomp = FALSE
} else {
  mcomp = TRUE
}

pcFiles <- c()
numPcFiles = length(chrvector)
for (a in 6:(6+numPcFiles)) {
  fullPath <- normalizePath(args[length(args)-a])
  pcFiles <- c(pcFiles, fullPath)
}

print(pcFiles)

diffcmp(chrvector, pcFiles, samplefile, lprmfile, res, mcomp)
#print("DiffCmp Run")
