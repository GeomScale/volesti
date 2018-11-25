V_poly_qhull <- function(P){
  
  if (P$t==2) {
    ps = convhulln(P$V,"FA")
    return(ps$vol)
  }
  
  ps1 = convhulln(P$V,"FA")
  ps2 = convhulln(P$V2,"FA")
  
  ps3 = convhulln(rbind(P$V,P$V2),"FA")
  
  volume = ps1$vol + ps2$vol - ps3$vol
  return(volume)
  
}