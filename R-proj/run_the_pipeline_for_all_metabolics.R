library(volesti)
library(Matrix)
library(Rmosek)

path_mets = "/home/tolis/data/metabolic_mat"
file.names <- dir(path_mets, pattern =".mat")

time_vec = c()

for (i in 10:length(file.names)) {
  
  path = paste0(path_mets, "/", file.names[i])
  print(path)
  
  tim = system.time({ res = apply_pipeline(path) })
  time_vec = c(time_vec, as.numeric(tim)[3])
  print(as.numeric(tim)[3])
  res$run_time = as.numeric(tim)[3]
  save(res, file = paste0(file.names[i],"_results.RData"))
  save(time_vec, file = "times_metabolics.RData")
  
}

save(time_vec, file = "times_metabolics.RData")
