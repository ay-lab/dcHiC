## This function helps to find any executable
.findExecutable <- function(exe, interactive=TRUE) {
  path <- Sys.which(exe)
  if(all(path=="")) {
    if(interactive) stop("Executable for ", paste(exe, collapse=" or "), " not found! Please make sure that the software is correctly installed and, if necessary, path variables are set.", call.=FALSE)
    return(character(0))
  }
  path[which(path!="")[1]]
}

cgse <- function(dcomp.file, gene_pos, organism='Homo.sapiens', go_enrich = FALSE, padj=1e-2, dzsc=3) {
  
  if (go_enrich==FALSE) {
    tmp <- tempfile(pattern = "temp", tmpdir = tempdir(), fileext = ".bed")
    dcomp<- read.table(dcomp.file, h=T, comment.char="", as.is=T)
    #dcomp<- dcomp[dcomp$padj < padj & dcomp$dZsc > dzsc,]
    colnames(dcomp)[1] <- "chr"
    dcomp <- dcomp[order(dcomp$chr),]
    rownames(dcomp) <- paste0(dcomp$chr,":",dcomp$start,"-",dcomp$end)
    write.table(dcomp[,1:3], file=tmp, col.names=F, row.names=F, sep="\t", quote=F)
    cmd <- paste0(.findExecutable("bedtools")," map -a ",tmp," -b ",gene_pos," -c 4 -o collapse -null NA")
    df <- read.table(text=system(cmd, wait=T,intern=T),h=F)
    rownames(df) <- paste0(df$V1,":",df$V2,"-",df$V3)
    dcomp[,"gene"] <- df[rownames(dcomp),4]
    write.table(na.omit(dcomp), file=paste0("Genes.",basename(dcomp.file)), row.names=F, sep="\t", quote=F)
    
    stop("Stopped for cGSEA. Not an error. Change 'go_enrich' field in cgsea.r if you would like to perform GO analysis instead.")
  
  }
  
  else {
    library(GOfuncR)
    tmp <- tempfile(pattern = "temp", tmpdir = tempdir(), fileext = ".bed")
    dcomp<- read.table(dcomp.file, h=T, comment.char="", as.is=T)
    dcomp<- dcomp[dcomp$padj < padj & dcomp$dZsc > dzsc,]
    colnames(dcomp)[1] <- "chr"
    dcomp <- dcomp[order(dcomp$chr),]
    rownames(dcomp) <- paste0(dcomp$chr,":",dcomp$start,"-",dcomp$end)
    write.table(dcomp[,1:3], file=tmp, col.names=F, row.names=F, sep="\t", quote=F)
    cmd <- paste0(.findExecutable("bedtools")," map -a ",tmp," -b ",gene_pos," -c 4 -o collapse -null NA")
    df <- read.table(text=system(cmd, wait=T,intern=T),h=F)
    rownames(df) <- paste0(df$V1,":",df$V2,"-",df$V3)
    dcomp[,"gene"] <- df[rownames(dcomp),4]
    write.table(na.omit(dcomp), file=paste0("Genes.",basename(dcomp.file)), row.names=F, sep="\t", quote=F)

    df <- na.omit(df)
    gene_names <- as.vector(df[,4])
    gene_names <- strsplit(paste0(unique(gene_names[!is.na(gene_names)]),collapse=","),",")[[1]]
    gene_names <- gene_names[gene_names != "NA"]
   
    input_hyper <- data.frame(gene_ids=gene_names, is_candidate=1)
    res_hyper <- go_enrich(input_hyper, n_randset=100, organismDb=organism)
    mf <- res_hyper[[1]][res_hyper[[1]]$ontology == "molecular_function" & res_hyper[[1]]$FWER_overrep <= 0.01,]
    bp <- res_hyper[[1]][res_hyper[[1]]$ontology == "biological_process" & res_hyper[[1]]$FWER_overrep <= 0.01,]
    rownames(mf) <- 1:nrow(mf)
    rownames(bp) <- 1:nrow(bp)
   
    data_mf <- list()
    for(i in 1:nrow(mf)) { 
      go_genes <- get_anno_genes(mf$node_id[i], database = organism)$gene
      common_gene <- intersect(go_genes, gene_names)
      regions <- list()
      k <- 1
      for(j in 1:nrow(df)) {
        if (length(grep(",",df$V4[j])) > 0) {
          region_genes <- strsplit(df$V4[j],",")[[1]] 
        } else {
          region_genes <- df$V4[j]
        }
        region_genes <- intersect(common_gene,region_genes)
        if (length(region_genes) > 0) {
          regions[[k]] <- rownames(df)[j] 
          k <- k + 1
        }
      }
      regions <- paste0(unique(unlist(regions)),collapse=",")
      print (head(regions))
      common_gene <- paste0(common_gene, collapse=",")
      data_mf[[i]] <- data.frame(mf[i,],genes=common_gene,regions)
    } 
    data_mf <- do.call(rbind, data_mf)
    write.table(data_mf, file="MF.txt", row.names=F, sep="\t", quote=F)
  
    data_bp <- list()
    for(i in 1:nrow(bp)) {
      go_genes <- get_anno_genes(bp$node_id[i], database = organism)$gene
      common_gene <- intersect(go_genes, gene_names)
      regions <- list()
      k <- 1
      for(j in 1:nrow(df)) {
        if (length(grep(",",df$V4[j])) > 0) {
          region_genes <- strsplit(df$V4[j],",")[[1]]
        } else {
          region_genes <- df$V4[j]
        }
        region_genes <- intersect(common_gene,region_genes)
        if (length(region_genes) > 0) {
          regions[[k]] <- rownames(df)[j]
          k <- k + 1
        }
      }
      regions <- paste0(unique(unlist(regions)),collapse=",")
      print (head(regions))
      common_gene <- paste0(common_gene, collapse=",")   
      data_bp[[i]] <- data.frame(bp[i,],genes=common_gene,regions)
    }
    data_bp <- do.call(rbind, data_bp)
    write.table(data_bp, file="BP.txt", row.names=F, sep="\t", quote=F)
  }
}

args = commandArgs(trailingOnly=TRUE) 
full_compt_details <- args[length(args)]
gene_bed <- args[length(args)-1]

cgse(full_compt_details, gene_bed)
#cgse("/mnt/BioAdHoc/Groups/vd-ay/abhijit/overflow/proj_overflow/Diffential_HiC_Compartment_Analysis/Test/Mouse_Bonev_Replicates/Bonev_data/DifferentialCompartment/NPC_vs_mESC_full_compartment_details.bedGraph", "mm10_gene_pos.bed", organism='Mus.musculus')
#cgse("NPC_vs_mESC_full_compartment_details.bedGraph", "mm10_gene_pos.bed")
