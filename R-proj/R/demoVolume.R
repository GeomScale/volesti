#' Run some experiments
#'
#' @return Print the computed volumes and the total time
#' @examples
#' testRvolEsti()
demoVolume <- function(){
  path=getwd()
  path=paste0(substr(path, start=1, stop=nchar(path)-7),'/data/')
  print(path)
  listofexamples=list.files(path)
  
  for(i in 1:length(listofexamples)){
    x=read.csv(paste0(path,listofexamples[i]))
    print(listofexamples[i])
    A=ineToMatrix(x)
    volume(list("matrix"=A,"test"=TRUE,"verbose"=TRUE))
  }
  
}

