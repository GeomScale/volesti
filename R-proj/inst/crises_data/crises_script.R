library(volesti)

starting_date = "2007-01-04"
stopping_date = "2010-01-04"

MatReturns <- read.table("https://stanford.edu/class/ee103/data/returns.txt", sep=",")

MatReturns = MatReturns[-c(1,2),]
dates = as.character(MatReturns$V1)
MatReturns = as.matrix(MatReturns[,-c(1,54)])
MatReturns = matrix(as.numeric(MatReturns [,]),nrow = dim(MatReturns )[1],
                    ncol = dim(MatReturns )[2], byrow = FALSE)

row1 = which(dates %in% starting_date)
row2 = which(dates %in% stopping_date)

indicators = compute_indicators(MatReturns =  MatReturns[row1:row2,], win_len = 60,
                                numSlices = 100, N = 1000000)
