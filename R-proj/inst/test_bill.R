library(volesti)

times_bil = numeric(length = 7)
counter = 1
for (w in c(1,5,seq(from=10,to=50, by=10))) {
  tim = system.time({ mat = control_billiard(200,200,200,w) })
  times_bil[counter] = as.numeric(tim)[3]
  counter = counter + 1
  print(w)
  save(times_bil, file = "times_bill200.RData")
}

save(times_bil, file = "times_bill200.RData")
