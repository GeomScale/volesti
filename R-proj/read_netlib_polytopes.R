name = '25fv47'  ##edit the name for another polytope
q=R.matlab::readMat(paste0('polytope_',name,'.mat'))

A = q$polytope[[1]]
b = q$polytope[[2]]
center = q$polytope[[3]]
radius = q$polytope[[4]]

