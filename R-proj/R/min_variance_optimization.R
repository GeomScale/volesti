  
  
  ############################################################################
  ### TOY EXAMPLE MINIMUM VARIANCE OPTMIZATION
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     02.11.2021
  # First version:    02.11.2021
  # --------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(quadprog)
  #wd <- ... # change working directory
  
  
  
  # --------------------------------------------------------------------------
  # LOAD DATA
  # --------------------------------------------------------------------------
  
  env <- readRDS("~/volume_approximation/R-proj/data/msci_ci.rds")
  X_est <- env$X_est
  
  
  
  
  # --------------------------------------------------------------------------
  # RUN OPTIMIZATION
  # --------------------------------------------------------------------------
  
  
  Dmat <- cov(X_est)
  
  Amat <- matrix(0, nrow = ncol(X_est) * 2 + 1, ncol = ncol(X_est),
                 dimnames = list(NULL, colnames(X_est)) )
  Amat[1, ] <- 1
  Amat[2:(ncol(X_est)+1), ] <- diag(ncol(X_est))
  Amat[(ncol(X_est)+2):nrow(Amat), ] <- -diag(ncol(X_est))
  Amat <- t(Amat)
  bvec <- c(1, rep(0, ncol(X_est)), rep(-1, ncol(X_est)) )
  dvec <- rep(0, ncol(X_est))
  
  opt <- quadprog::solve.QP( Dmat = Dmat,
                             dvec = dvec,
                             Amat = Amat,
                             bvec = bvec,
                             meq  = 1 )
  wghts <- setNames( opt$solution, colnames(X_est) )
  wghts
  
  barplot( wghts )
  
  
  
  
  