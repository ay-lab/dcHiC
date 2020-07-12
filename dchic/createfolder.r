createFolder <- function(name, changepath=0) {

  if (!file.exists(name)) {
    dir.create(name)
    cat (name," folder created\n")
  } else {
    cat (name," folder already exists\n")
  }
  if (changepath == 1) {
    setwd(name)
  }
}
