library(Rcpp)
library(hashmap)
source("hicmat.r")
source("createfolder.r")
sourceCpp("ijkmat.cpp")

args = commandArgs(trailingOnly=TRUE)
proj_name <- args[length(args)]
text_file <- args[length(args)-1]


df <- read.table(text_file, header=FALSE)
#matrices <- vector(mode="character", length=length(df[1]))
#bed <- vector(mode="character", length=length(df[2]))
#prefixes <- vector(mode="character", length=length(df[3]))

matrices <- c()
bed <- c()
prefixes <- c()

print(df)
# for (a in 1:length(df[1]))
#  matrices[a] <- df[[1]][[a]]
# for (a in 1:length(df[2]))
#  bed[a] <- df[[2]][[a]]
# for (a in 1:length(df[3]))
#  prefixes[a] <- df[[3]][[a]]

for (a in 1:nrow(df)) {
  loc <- as.character(df$V1[[a]])
  matrices <- c(matrices, loc)
}
for (a in 1:nrow(df)) {
  loc <- as.character(df$V2[[a]])
  bed <- c(bed, loc)
}
for (a in 1:nrow(df)) {
  loc <- as.character(df$V3[[a]])
  prefixes <- c(prefixes, loc)
}

class(matrices)
mat <- normalizePath(matrices)
bedf <- normalizePath(bed)

project <- proj_name
createFolder(project, 1)
hicmat(mat, bedf, prefixes)