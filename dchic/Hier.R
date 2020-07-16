library(FactoMineR)
args = commandArgs(trailingOnly=TRUE) 
path = args[length(args)]
source(file.path(path,"armcorrection.R"))
source(file.path(path,"selectpc.R"))

dataPath <- args[length(args)-1]
selectedGenome <- args[length(args)-2]
selectedChr <- args[length(args)-3]
selectedRes <- as.integer(args[length(args)-4])
selectedMode <- args[length(args)-5]
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
res.hmfa <- HMFA(comb_matrix, H = hierar, type = rep("c",length(hierar[[1]])), ncp = 2, name.group = names, graph=FALSE)
pcNum <- selectpc(res.hmfa, selectedChr)
sink("pc_decision.txt")
cat(pcNum)
sink()
if (selectedMode != "None") {
  if (dataPath == "None") {
    res.hmfa <- armcorrection(obj=res.hmfa, genome=selectedGenome, pc=pcNum, chr=selectedChr, resolution=selectedRes, choice=selectedMode)
  } else {
    res.hmfa <- armcorrection(obj=res.hmfa, genome=selectedGenome, pc=pcNum, chr=selectedChr, resolution=selectedRes, choice=selectedMode, folder = normalizePath(dataPath))
  }
}

print(currentDir)
setwd(currentDir)

saveRDS(res.hmfa, "Rsession.rds")

bins <- labels(res.hmfa[["partial"]][[1]])

bigframe <- vector("list", numElements)
pc2 <- vector("list", numElements)

for (i in 1:(numElements*2)) {
  if ((i+1)%%2 == 0) {
    startpos = (i-1)*size+1
    endpos = startpos+size-1 # change 2/16; never had prior effects b/c this is odd
    vec <- c()
    for (j in startpos:endpos) {
      vec <- c(vec, res.hmfa[["partial"]][[1]][[j]])
    }
    bigframe[[(i+1)/2]] <- vec
  }
}


# PC 2 Stuff
for (i in 1:(numElements*2)) {
  if ((i)%%2 == 0) {
    startpos = (i-1)*size+1
    endpos = startpos+size-1 #change 2-16
    vec <- c()
    #print(startpos)
    #print(endpos)
    for (j in startpos:endpos) {
      vec <- c(vec, res.hmfa[["partial"]][[1]][[j]])
    }
    pc2[[i/2]] <- vec
  }
}

for (i in 1:numElements) {
  str_init <- paste("hmfa_chrRAW_", general_names[[i]], sep = "") #numChr is a legacy feature
  str_init1 <- paste(str_init, "_exp_", sep = "")
  str_res <- paste(str_init1, i, sep = "")
  str_res2 <- paste(str_res, ".txt", sep= "")
  exp_name <- paste("exp", i, sep = "")
  

  file.create(str_res2)
  fileConn<-file(str_res2)
  vec <- c()
  vec <- c(vec, paste("\t", "\t", "Dim.1"))
  for (j in 1:size) {
    vec <- c(vec, paste(bins[[1]][[j]], "\t", bigframe[[i]][[j]]))
    writeLines(vec, fileConn)
  }
  close(fileConn)
}

for (i in 1:numElements) {
  str_init <- paste("pc2_chrRAW_", general_names[[i]], sep = "") #numChr is a legacy feature
  str_init1 <- paste(str_init, "_exp_", sep = "")
  str_res <- paste(str_init1, i, sep = "")
  str_res2 <- paste(str_res, ".txt", sep= "")
  exp_name <- paste("exp", i, sep = "")
  
  file.create(str_res2)
  fileConn<-file(str_res2)
  vec <- c()
  vec <- c(vec, paste("\t", "\t", "Dim.2"))
  for (j in 1:size) {
    vec <- c(vec, paste(bins[[1]][[j]], "\t", pc2[[i]][[j]]))
    writeLines(vec, fileConn)
  }
  close(fileConn)
}

sink("coordinates.txt")
for (i in 1:2) { #this is assuming you only ever have two levels
  print(res.hmfa[["group"]][["coord"]][[i]])
}
sink()

