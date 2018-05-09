# dataReading accepts a datapath for an uploaded file,
# and outputs a data matrix based on the upload contents.
# Functionality with Tab-separated files removed.

dataReading <- function(upload.path) {
  if(file.exists("new.txt")) {file.remove("new.txt")}
  write("mass,intensity", file = "new.txt", append = TRUE)
  con <- file(upload.path, open = 'r')
  while(TRUE) {
    line <- readLines(con, n = 1)
    # lines starting with hashes or containing 0 values for intensity are removed.
    if(length(line) == 0) break
    else if(startsWith(line, "#") | endsWith(line, ",0") | endsWith(line, ",0.00") | endsWith(line, ",0.0000")) {
    } else {
      write(line, file = "new.txt", append = TRUE)
    }
  }
  
  data <- read.table("new.txt",header=TRUE,sep = ",")
  
  if(ncol(data) != 2) {
    data <- matrix(c(0,0,0,0), ncol=2, nrow=2)
  }
  
  file.remove("new.txt")
  msMatrix <- as.matrix(data)
  return(msMatrix)
}