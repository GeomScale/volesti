if(FALSE) {
path_mets = "/home/tolis/data/metabolic"
file.names <- dir(path_mets, pattern =".mat")

for (i in 1:length(file.names)) {
  
  path = paste0(path_mets, "/", file.names[i])
  print(path)
  
  res = apply_pipeline(path)
  print(res)
  
}
}
