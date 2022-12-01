#!/usr/bin/env Rscript

library(optparse)

## The following code is copied from Tcseq package
## Wu M, Gu L (2021). TCseq: Time course sequencing data analysis. R package version 1.18.0.
## Modified as per requirement by Abhijit Chakraborty

clust <- methods::setClass("clust", slots = c(method = "character",
                                     dist = "character",
                                     data = "matrix",
                                     centers = "matrix",
                                     cluster = "integer",
                                     membership = "matrix"))

pcaclust <- function(x, algo, k, dist = "euclidean", centers = NULL,
                      standardize = TRUE, ...) {
  if (is.matrix(x)) {
    data.tmp <- x
  }else{
    data.tmp <- x@tcTable
  }
  if (standardize) {
    for (i in seq_len(nrow(data.tmp))) {
      data.tmp[i, ] <- (data.tmp[i, ] - mean(data.tmp[i, ], na.rm = TRUE))/sd(data.tmp[i, ], na.rm = TRUE)
    }
    data.tmp <- data.tmp[complete.cases(data.tmp), ]
  }
  object <- methods::new("clust")
  object@method <- algo
  object@dist <- dist
  object@data <- data.tmp
  
  res <- .pcaclust(data = data.tmp, algo = algo, k = k, dist = dist,
                    centers = centers, ...)
  
  if (algo == "cm") {
    object@cluster <- res$cluster
    object@membership <- res$membership
    object@centers <- res$centers
  } else {
    object@cluster <- res$cluster
    object@centers <- res$centers
  }
  if (is.matrix(x)) {
    object
  } else {
    x@clusterRes <- object
    x
  }
}

# perform time course clustering
.pcaclust <- function(data, algo, k, centers = NULL,
                       dist = "euclidean", ...) {
  if (!algo %in% c("pam", "km", "hc", "cm")) {
    stop("clustering method should be one of 'pam','km','hc','cm'")
  }
  if (!dist %in% c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) {
    stop("Distance metric should be 'correlation', or one of the distance measures in dist function")
  }
  if (algo == "km") {
    if(dist != "euclidean"){
      stop("kmeans only support euclidean metric; for other distance metrices, please see the help page")
    }
  }
  if (algo == "cm" ) {
    if(!dist %in% c("euclidean", "manhattan")){
      stop("cmeans only support euclidean or mahattan distance metrics")
    }
  }
  
  d <- NULL
  if (algo %in% c("pam", "hc")) {
    if (dist == "correlation") {
      d <- as.dist(1 - cor(t(data)))
    }
    if (dist != "correlation") {
      d <- dist(data, method = dist)
    }
  }
  clustres <- list()
  if (algo != "hc") {
    if (!is.null(centers)) {
      if (nrow(centers) != k) {
        stop("Number of rows of centers must be equal to k")
      }
    }
  }
  clustres <- switch(algo, km = {
    if (!is.null(centers)) {
      res <- kmeans(data, centers = centers, ...)
    } else {
      res <- kmeans(data, centers = k, ...)
    }
    clustres$cluster <- res$cluster
    clustres$centers <- res$centers
    clustres
  }, pam = {
    if (!is.null(centers)) {
      ind <- data[, 1] %in% centers[, 1]
      ind <- which(ind)
      if (length(ind) != k) {
        stop("For 'pam', centers must be chosen from the data")
      } else {
        res <- pam(d, k = k, medoids = ind, ...)
      }
    }
    res <- pam(d, k = k, ...)
    clustres$cluster <- res$clustering
    clustres$centers <- data[res$medoids, ]
    clustres
  }, hc = {
    tree <- hclust(d, ...)
    res <- cutree(tree, k = k)
    clustres$cluster <- res
    clustres$centers <- matrix(0, 0, 0)
    clustres
  }, cm = {
    if (!is.null(centers)) {
      res <- e1071::cmeans(data, centers = centers, ...)
    } else {
      res <- e1071::cmeans(data, centers = k, ...)
    }
    clustres$cluster <- res$cluster
    clustres$centers <- res$centers
    clustres$membership <- res$membership
    clustres
  })
  clustres
}

pcaclustplot <- function(object = NULL, categories = "timepoint",
                          value = "expression", cols = NULL,
                          cl.color = "gray50",
                          membership.color = rainbow(30, s = 3/4, v = 1, start = 1/6),
                          title.size = 18, axis.line.size = 0.6,
                          axis.title.size = 18,
                          axis.text.size = 16, legend.title.size = 14,
                          legend.text.size = 14, outprefix= "pca_clusters") {

  if (class(object) != "clust" && class(object) != "TCA") {
    stop("object should be a 'timeclust' object or a 'TCA' object")
  }
  if (class(object) == "clust") {
    data <- object@data
    cluster <- object@cluster
    membership <- object@membership
  }
  if (class(object) == "TCA") {
    data <- object@clusterRes@data
    cluster <- object@clusterRes@cluster
    membership <- object@clusterRes@membership
  }
  ncl <- max(cluster)
  membercolor <- vector(length = length(cluster))
  membervalue <- list()
  counter <- 0
  if (!sum(dim(membership) == 0) == 2) {
    color <- membership.color
    colorseq <- seq(0, 1, length = length(color))
    for (i in seq_len(ncl)) {
      mtmp <- membership[cluster == i, i]
      membervalue[[i]] <- mtmp
      for (j in seq_len(length(mtmp))) {
        counter <- counter + 1
        ind <- which(abs(colorseq - mtmp[j]) == min(abs(colorseq - mtmp[j])))
        membercolor[counter] <- color[ind]
      }
    }
    membervalue <- unlist(membervalue)
    names(membercolor) <- membervalue
  }

  plotlist <- list()
  for (i in seq_len(ncl)) {
    title <- paste0("Cluster ", i)
    dtmp <- data[cluster == i, ]
    a <- which(cluster == i)
    if (length(a) == 1) {
      dtmp <- data.frame(time = 1:length(dtmp), value = dtmp)
      if (!sum(dim(membership) == 0) == 2) {
        m <- membership[cluster == i, i]
        colorname = toString(m)
        plotlist[[i]] <- ggplot2::ggplot(dtmp, ggplot2::aes(x = time, y = value)) +
          ggplot2::geom_line(colour = membercolor[colorname]) + ggplot2::theme_bw() +
          ggplot2::ggtitle(title) +
          ggplot2::scale_x_continuous(breaks = dtmp$time,
                             labels = row.names(dtmp)) +
          ggplot2::labs(x = categories, y = value) +
          ggplot2::theme(plot.title = ggplot2::element_text(size = title.size),
                axis.line.x = ggplot2::element_line(color = "black",
                                           size = axis.line.size),
                axis.line.y = ggplot2::element_line(color = "black",
                                           size = axis.line.size),
                axis.title = ggplot2::element_text(size = axis.title.size),
                axis.text = ggplot2::element_text(size = axis.text.size),
                legend.position = "none", panel.border = ggplot2::element_blank(),
                panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank())
      } else {
        plotlist[[i]] <- ggplot2::ggplot(dtmp, ggplot2::aes(x = time, y = value)) +
          ggplot2::geom_line(colour = cl.color) + ggplot2::theme_bw() + ggplot2::ggtitle(title) +
          ggplot2::scale_x_continuous(breaks = dtmp$time,
                             labels = row.names(dtmp)) +
          ggplot2::labs(x = categories, y = value) +
          ggplot2::theme(plot.title = ggplot2::element_text(size = title.size),
                axis.line.x = ggplot2::element_line(color = "black",
                                           size = axis.line.size),
                axis.line.y = ggplot2::element_line(color = "black",
                                           size = axis.line.size),
                axis.title = ggplot2::element_text(size = axis.title.size),
                axis.text = ggplot2::element_text(size = axis.text.size),
                legend.position = "none", panel.border = ggplot2::element_blank(),
                panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank())
      }
    } else {
      dtmp_m <- reshape2::melt(dtmp)
      colnames(dtmp_m) <- c("group", "time", "value")
      if (sum(dim(membership) == 0) == 2) {
        plotlist[[i]] <- ggplot2::ggplot(dtmp_m, ggplot2::aes(x = time, y = value)) +
          ggplot2::geom_line(ggplot2::aes(group = group), colour = cl.color) +
          ggplot2::theme_bw() + ggplot2::ggtitle(title) +
          ggplot2::labs(x = categories, y = value) +
          ggplot2::theme(plot.title = ggplot2::element_text(size = title.size),
                axis.line.x = ggplot2::element_line(color = "black",
                                           size = axis.line.size),
                axis.line.y = ggplot2::element_line(color = "black",
                                           size = axis.line.size),
                axis.title = ggplot2::element_text(size = axis.title.size),
                axis.text = ggplot2::element_text(size = axis.text.size),
                legend.position = "none", panel.border = ggplot2::element_blank(),
                panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank())
      }
      if (!sum(dim(membership) == 0) == 2) {
        mem <- membership[cluster == i, i]
        mem1 <- data.frame(group = names(mem), member = mem)
        dtmp_m1 <- merge(dtmp_m, mem1, by = "group")
        colnames(dtmp_m1) <- c("group", "time", "value", "membership")
        dtmp_m1 <- dtmp_m1[order(dtmp_m1[, 4]), ]
        new.factor <- unique(as.vector(dtmp_m1$group))
        dtmp_m1$group <- factor(dtmp_m1$group, levels = new.factor)

        plotlist[[i]] <- ggplot2::ggplot(dtmp_m1, ggplot2::aes(x = time, y = value,
                                             colour = membership)) +
          ggplot2::geom_line(ggplot2::aes(group = group)) +
          ggplot2::scale_colour_gradientn(colours = membership.color) +
          ggplot2::guides(colour = ggplot2::guide_colourbar()) + ggplot2::theme_bw() +
          ggplot2::ggtitle(title) + ggplot2::labs(x = categories, y = value) +
          ggplot2::theme(plot.title = ggplot2::element_text(size = title.size),
                axis.line.x = ggplot2::element_line(color = "black",
                                           size = axis.line.size),
                axis.line.y = ggplot2::element_line(color = "black",
                                           size = axis.line.size),
                axis.title = ggplot2::element_text(size = axis.title.size),
                axis.text = ggplot2::element_text(size = axis.text.size),
                legend.title = ggplot2::element_text(size = legend.title.size),
                legend.text = ggplot2::element_text(size = legend.title.size),
                panel.border = ggplot2::element_blank(),
                panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank())


      }
    }
  }
  plots <- c(list, plotlist)
  for(i in 1:length(plots)) {
    print(plots[[i]])
  }
  dev.off() 
}

#Get the options
option_list = list(
  make_option(c("--pcafile"), type="character", help="dcHiC output file for example the differential.<intra/inter>_sample_group.Filtered.pcQnm.bedGraph file\n"),
  make_option(c("--samplefile"), type="character", help="A single column file that defines which columns to select from the bedGraph file\n"),
  make_option(c("--kcenter"), type="integer", default=6, help="The number of centers to find from bedGraph file [Default 6]\n"),
  make_option(c("--genebed"), type="character", help="A sorted bed file with gene names in its fourth column\n"),
  make_option(c("--minprob"), type="numeric", default=0, help="Minimum probability for cluster membership [Default 0]\n"),
  make_option(c("--output"), type="character", default="pc_clustered", help="The prefix for all the output files [Default pc_clustered]\n\n
  \t##### Note #####
  \tPlease cite Wu M, Gu L (2021). TCseq: Time course sequencing data analysis. R package version 1.18.0 if you're using this code along with dcHiC\n")
)

opt <- parse_args(OptionParser(option_list=option_list))

pcafile    <- as.character(opt$pcafile)
samplefile <- as.character(opt$samplefile)
kcenter    <- as.integer(opt$kcenter)
genebed    <- as.character(opt$genebed)
minprob    <- as.numeric(opt$minprob)
output     <- as.character(opt$output)

pcafile <- read.table(pcafile, h=T, as.is=T)
samplefile <- read.table(samplefile, h=F, as.is=T)

pcafile <- pcafile[,c("chr","start","end",as.character(samplefile$V1))]
pcafile <- pcafile[order(pcafile$chr, pcafile$start),]
rownames(pcafile) <- paste0(pcafile$chr,"_",pcafile$start,"_",pcafile$end)
pcafile <- as.matrix(pcafile[,-c(1:3)])
head(pcafile)

pcacl <- pcaclust(pcafile, algo="cm", k=kcenter, standardize = TRUE)
pcaclustplot(pcacl, value="PCA (Quantile Normalized & Z-score)", cols=1)
system(paste0("mv Rplots.pdf Cluster_plots.",output,".pdf"), wait=T)

cl_df <- as.data.frame(pcacl@membership)
cl_df <- cl_df[apply(cl_df, 1, max) > minprob,]
cat("Keeing members that are above minimum probability of ",minprob,"\n")
cl_df[,"cluster"] <- paste0("Cluster_",apply(cl_df, 1, which.max))

write.table(pcacl@centers, file=paste0("Centers.",output,".txt"), sep="\t", quote=F)
write.table(cl_df, file=paste0("Membership.",output,".txt"), sep="\t", quote=F)
cl_df[,c("chr","start","end")] <- data.frame(do.call(rbind,strsplit(rownames(cl_df),"_")))
head(cl_df)

for(i in 1:kcenter) {
  write.table(cl_df[cl_df$cluster==paste0("Cluster_",i),c("chr","start","end","cluster")], file=paste0("Cluster_",i,".",output,".bedGraph"), row.names=F, col.names=F, sep="\t", quote=F)
  bedmap <- data.table::fread(paste0("bedtools map -a Cluster_",i,".",output,".bedGraph -b ",genebed," -c 4 -o distinct -null NA"))
  bedmap <- na.omit(bedmap)
  genes  <- sort(unique(unlist(strsplit(bedmap$V5,","))))
  write.table(genes, file=paste0("Genes_Cluster_",i,".",output,".txt"), row.names=F, col.names=F, sep="\t", quote=F)

}
