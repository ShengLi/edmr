# Read text file fast by predicting and feeding \code{colClass} object
readTableFast<-function(filename,header=T,skip=0,sep="", ...)
{
  tab5rows <- read.table(filename, header = header,skip=skip,sep=sep, nrows = 100, ...)
  classes  <- sapply(tab5rows, class)
  return( read.table(filename, header = header,skip=skip,sep=sep, colClasses = classes, ...)  )
}

# split function for strings
splitn=function (strings, field, n)
{
  sapply(strsplit(strings, field), "[[", n)
}
