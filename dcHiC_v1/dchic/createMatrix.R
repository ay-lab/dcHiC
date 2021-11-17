library(Rcpp)
library(hashmap)
library(base)

thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}

args = commandArgs(trailingOnly=TRUE)
source(file.path(dirname(thisFile()),"hicmat.r"))
source(file.path(dirname(thisFile()),"createfolder.r"))
sourceCpp(file.path(dirname(thisFile()),"ijkmat.cpp"))


proj_name <- args[length(args)]
text_file <- args[length(args)-1]

df <- read.table(text_file, header=FALSE)

matrices <- c()
bed <- c()
prefixes <- c()

print(df)

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
