library(volesti)

load("~/volume_approximation/R-proj/inst/opt_sols.RData")

dims=seq(from=30,to=60,by=10)
times_hmc = numeric(length = 10)
iters_hmc = numeric(length = 10)

times_hmc = matrix(0,length(dims),1)
iters_hmc = matrix(0,length(dims),1)

times_hnr = matrix(0,length(dims),1)
iters_hnr = matrix(0,length(dims),1)

dim_it=0
for (d in dims) {
  dim_it = dim_it + 1

  for (i in 1:1) {
    
    iters = 0
    avg_time = 0
    #for (j in 1:1) {
      
#      tim = system.time({ iters = iters + sdp_approx(d, i+1, opt_sols[dim_it, i+1], 0.05, 1, 1) })
#      avg_time = avg_time + as.numeric(tim)[3]
    
    #}
    #iters_hmc[dim_it, i] = iters / 1
    #times_hmc[dim_it, i] = avg_time / 1
    
    
    iters = 0
    avg_time = 0
    for (j in 1:1) {
      
      
      tim = system.time({ iters = iters + sdp_approx(d, i+1, opt_sols[dim_it, i+1],0.05, 2, 4*sqrt(d)) })
      avg_time = avg_time + as.numeric(tim)[3]
      
    }
    iters_hnr[dim_it, i] = iters / 1
    times_hnr[dim_it, i] = avg_time / 1
    
    save(iters_hmc, times_hmc, iters_hnr, times_hnr, file = "sdps_iters_times_small01.RData")
  
  }
}