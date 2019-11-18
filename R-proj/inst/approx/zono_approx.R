library(volesti)

d=10
orders = c(5,10,15,20,50)

times10=c()
vols10=c()
vols_red10=c()
ratios10=c()
for (i in orders) {
  P = GenRandZonotope(d, d*i, dist = "uniform")
  tim = system.time({ vec = zono_approx(P, random_walk = "BilW", fit_ratio = TRUE) })
  
  times10 = c(times10, as.numeric(tim[3]))
  vols10 = c(vols10, vec$vol)
  vols_red10 = c(vols_red10, vec$vol_red)
  ratios10 = c(ratios10, vec$fit_ratio)
}

save(times10, vols10, vols_red10, ratios10, file = "d10_zono_approx.RData")


d=15
times15=c()
vols15=c()
vols_red15=c()
ratios15=c()
for (i in orders) {
  P = GenRandZonotope(d, d*i, dist = "uniform")
  tim = system.time({ vec = zono_approx(P, random_walk = "BilW", fit_ratio = TRUE) })
  
  times15 = c(times15, as.numeric(tim[3]))
  vols15 = c(vols15, vec$vol)
  vols_red15 = c(vols_red15, vec$vol_red)
  ratios15 = c(ratios15, vec$fit_ratio)
}

save(times15, vols15, vols_red15, ratios15, file = "d15_zono_approx.RData")



d=20
times20=c()
vols20=c()
vols_red20=c()
ratios20=c()

P = GenRandZonotope(d, 200, dist = "uniform")
tim = system.time({ vec = zono_approx(P, random_walk = "BilW", fit_ratio = TRUE) })
  
times20 = c(times20, as.numeric(tim[3]))
vols20 = c(vols20, vec$vol)
vols_red20 = c(vols_red20, vec$vol_red)
ratios20 = c(ratios20, vec$fit_ratio)
  
P = GenRandZonotope(d, 300, dist = "uniform")
tim = system.time({ vec = zono_approx(P, random_walk = "BilW", fit_ratio = TRUE) })
  
times20 = c(times20, as.numeric(tim[3]))
vols20 = c(vols20, vec$vol)
vols_red20 = c(vols_red20, vec$vol_red)
ratios20 = c(ratios20, vec$fit_ratio)

save(times20, vols20, vols_red20, ratios20, file = "d20_zono_approx.RData")



d=30
times30=c()
vols30=c()
vols_red30=c()
ratios30=c()

P = GenRandZonotope(d, 300, dist = "uniform")
tim = system.time({ vec = zono_approx(P, random_walk = "BilW", fit_ratio = TRUE) })

times30 = c(times30, as.numeric(tim[3]))
vols30 = c(vols30, vec$vol)
vols_red30 = c(vols_red30, vec$vol_red)
ratios30 = c(ratios30, vec$fit_ratio)

P = GenRandZonotope(d, 450, dist = "uniform")
tim = system.time({ vec = zono_approx(P, random_walk = "BilW", fit_ratio = TRUE) })

times30 = c(times30, as.numeric(tim[3]))
vols30 = c(vols30, vec$vol)
vols_red30 = c(vols_red30, vec$vol_red)
ratios30 = c(ratios30, vec$fit_ratio)

save(times30, vols30, vols_red30, ratios30, file = "d30_zono_approx.RData")

