library(volesti)


N = 10
dims = seq(from = 20, to = 200, by=20)

vol_vec=matrix(0,length(dims),N)
steps_vec = matrix(0,length(dims),N)
time_vec=matrix(0,length(dims),N)
count = 1
for (d in dims) {
  print(d)
  generator_sdp(d,d)
  
  for (i in 1:N) {
    print(i)
    tim = system.time({ vol = volume(paste0('sdp_prob_',d,'_',d,'.txt')) })
    
    vol_vec[count,i] = vol[1]
    steps_vec[count,i] = vol[2]
    time_vec[count,i] = as.numeric(tim)[3]
    
    save(vol_vec, file = "vol_cb_spectra_40_200.RData")
    save(time_vec, file = "times_cb_spectra_40_200.RData")
    save(steps_vec, file = "steps_cb_spectra_40_200.RData")
  }
  count = count + 1
}