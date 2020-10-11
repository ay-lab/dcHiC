library(FactoMineR)
args = commandArgs(trailingOnly=TRUE) 
path = args[length(args)]
source(file.path(path,"pcselect.r"))

ncp_user <- as.numeric(args[length(args)-1])
dataPath <- args[length(args)-2]
selectedGenome <- args[length(args)-3]
selectedChr <- args[length(args)-4]
selectedRes <- as.integer(args[length(args)-5])
size = as.numeric(args[length(args)-6])
numElements = as.numeric(args[length(args) - 7])
numGroups = as.numeric(args[length(args) - 8])

comb_matrix <- read.table(args[1])
for (i in 2:numElements) {
  comb_matrix <- cbind(comb_matrix, read.table(args[i]))
}

hierar <- vector("list", 2)
hierar[[1]] <- rep(size, numElements)

pos = numElements + 1
group_vector <- c()
for (i in pos:(pos+numGroups-1)) {
  group_vector <- c(group_vector,as.numeric(args[i]))
}

hierar[[2]] <- group_vector

print("Hierarchy")
print(hierar)

names <- vector("list", 2)

general_names <- vector("list", numElements)
for (i in (pos+numGroups):(pos+numGroups+numElements-1)) {
  general_names[[i-pos-numGroups+1]] <- args[i]
}

group_names <- vector("list", numGroups)
for (i in (pos+numGroups+numElements):(pos+numGroups+numElements+numGroups-1)) {
  group_names[[i-pos-numGroups-numElements+1]] <- args[i] 
}

names[[1]] <- general_names
names[[2]] <- group_names

currentDir <- getwd()
res.hmfa <- HMFA(comb_matrix, H = hierar, type = rep("c",length(hierar[[1]])), ncp = ncp_user, name.group = names, graph=FALSE)

name <- paste("Rsession_HMFAObj_", selectedChr, ".rds", sep = "")
saveRDS(res.hmfa, name)

if (dataPath == "None") {
  pc.res <- pcselect(obj=res.hmfa, chr=selectedChr, genome=selectedGenome, resolution=selectedRes)
} else {
  pc.res <- pcselect(obj=res.hmfa, chr=selectedChr, genome=selectedGenome, resolution=selectedRes, folder = normalizePath(dataPath))
}

print(currentDir)
setwd(currentDir)

name <- paste("Rsession_", selectedChr, ".rds", sep = "")
saveRDS(res.hmfa, name)

bins <- labels(res.hmfa[["partial"]][[1]])

for (i in 1:numElements) {
  col <- c("\t\tDim.1")
  str_init <- paste("hmfa_chrRAW_", general_names[[i]], sep = "") #numChr is a legacy feature
  str_init1 <- paste(str_init, "_exp_", sep = "")
  str_res <- paste(str_init1, i, sep = "")
  str_res2 <- paste(str_res, ".txt", sep= "")
  exp_name <- paste("exp", i, sep = "")
  PCvals <- as.data.frame(pc.res[[i]]$pc, row.names = bins[[1]])
  write.table(PCvals, file = str_res2, row.names = T, quote = F, col.names = col)
}

sink("coordinates.txt")
for (i in 1:2) { # this is assuming you only ever have two levels
  print(res.hmfa[["group"]][["coord"]][[i]])
}
sink()
