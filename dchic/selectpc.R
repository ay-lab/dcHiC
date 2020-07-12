selectpc <- function(hmfa.obj, chr=NA, span=0.25) {

  pc1 <- hmfa.obj$ind[[1]][,1]
  pc2 <- hmfa.obj$ind[[1]][,2]

  #pc1_val <- pc1
  #pc2_val <- pc2 

  pc1[pc1 > 0] <- 1
  pc1[pc1 <= 0] <- -1
  pc2[pc2 > 0] <- 1
  pc2[pc2 <= 0] <- -1
  df <- data.frame(bin=c(1:length(pc1)), pc1, pc2)
  pc1_id <- 1
  pc2_id <- 1
  df[,"pc1_id"] <- 0
  df[,"pc2_id"] <- 0
  df$pc1_id[1]  <- pc1_id
  df$pc2_id[1]  <- pc2_id
  i  <- 2
  while (i <= nrow(df)) {
    pc1_sign <- df$pc1[(i-1)]
    pc2_sign <- df$pc2[(i-1)]
    if (df$pc1[i] == pc1_sign) {
      df$pc1_id[i] <- pc1_id
    } else if (df$pc1[i] != pc1_sign) {
      pc1_id <- pc1_id + 1
      df$pc1_id[i] <- pc1_id
    }

    if (df$pc2[i] == pc2_sign) {
      df$pc2_id[i] <- pc2_id
    } else if (df$pc2[i] != pc2_sign) {
      pc2_id <- pc2_id + 1
      df$pc2_id[i] <- pc2_id
    }
    i <- i + 1
  }

  df_pc1_min_coord <- aggregate(bin ~ pc1_id, min, data=df)$bin
  df_pc1_max_coord <- aggregate(bin ~ pc1_id, max, data=df)$bin
  df_pc1_sum <- aggregate(pc1 ~ pc1_id, sum, data=df)$pc1

  df_pc2_min_coord <- aggregate(bin ~ pc2_id, min, data=df)$bin
  df_pc2_max_coord <- aggregate(bin ~ pc2_id, max, data=df)$bin  
  df_pc2_sum <- aggregate(pc2 ~ pc2_id, sum, data=df)$pc2

  df_pc1 <- data.frame(df_pc1_min_coord,df_pc1_max_coord,df_pc1_sum)
  df_pc2 <- data.frame(df_pc2_min_coord,df_pc2_max_coord,df_pc2_sum)
 
  df_pc1 <- df_pc1[abs(df_pc1$df_pc1_sum) > (nrow(df) * span),]
  df_pc2 <- df_pc2[abs(df_pc2$df_pc2_sum) > (nrow(df) * span),]
  
  if (nrow(df_pc1) == nrow(df_pc2)) {
    pc <- 1
  } else if (nrow(df_pc1) < nrow(df_pc2)) {
    pc <- 1
  } else if (nrow(df_pc1) > nrow(df_pc2)) {
    pc <- 2
  }
  
  return(pc)
  
  #plot(pc1_val, type="l", col="red", main=paste0(chr," PC selected : ",pc), ylim=c(-10,10))
  #lines(pc2_val, col="blue")
  #abline(h=0)
  #legend(10, 8, legend=c("PC 1", "PC 2"), col=c("red", "blue"), lty=1, cex=0.8)
}

# args = commandArgs(trailingOnly=TRUE) 
# chrNum <- as.integer(args[length(args)])
# obj <- readRDS("Rsession.rds")
# 
# result <- selectpc(obj, chrNum)
# sink("pc_decision.txt")
# cat(result)
# sink()
# cat("Chromosome ", paste0(chrNum), " Decision: PC ", result)

#pdf("PC_selection.pdf")
#chr <- paste0("chr",c(1:22,"X"))
#for(i in 1:length(chr)) {
#  obj <- readRDS(paste0("Replicate_Analysis/RDSdata/",chr[i],".hmfaData.RDS"))
#  cat (paste0(chr[i]),"\t",selectpc(obj, chr[i]),"\n")
#}
#dev.off()
