library(volesti)
library(Matrix)
library(Rmosek)

path_mets = "/home/tolis/data/metabolic_mat/mat_socg/homo"
file.names <- dir(path_mets, pattern =".mat")

time_vec = c()

for (i in 1:length(file.names)) {
  
  path = paste0(path_mets, "/", file.names[i])
  print(path)
  
  tim = system.time({ res = apply_pipeline(path, remove_biomass = FALSE, save_files = TRUE) })
  time_vec = c(time_vec, as.numeric(tim)[3])
  print(as.numeric(tim)[3])
  res$run_time = as.numeric(tim)[3]
  save(res, file = paste0(file.names[i],"_results.RData"))
  save(time_vec, file = "times_metabolics.RData")
  
}

save(time_vec, file = "times_metabolics.RData")
