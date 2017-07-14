

read.data <- function(fname) {
  data = read.csv(fname, header=FALSE)
  rnames = as.matrix(data[-1,1])
  cnames = as.matrix(data[1,-1])
  data = data[-1,]
  data = as.matrix(data[,-1])
  data = matrix(as.numeric(data), dim(data))
  rownames(data) = rnames
  colnames(data) = cnames

  return(data)
}
