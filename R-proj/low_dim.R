# calculate the Moore-Penrose pseudo-inverse of matrix
pseudoinverse <- function (m) {
  msvd <- svd(m)
  if (length(msvd$d) == 0) {
    array(0, dim(m)[2:1])
  } else {
    s <- matrix(0, nrow=length(msvd$d), ncol=length(msvd$d))
    diag(s)[msvd$d != 0] <- 1/msvd$d[msvd$d != 0]
    msvd$v %*% (s %*% t(msvd$u))
  }
}



# generate basis for solution space of equality constraints
solution.basis <- function(constr) {
  stopifnot(all(constr$dir == '='))
  A <- constr$constr
  b <- constr$rhs
  n <- ncol(A)
  
  A.inv <- pseudoinverse(A)
  
  if (isTRUE(all.equal(diag(n), A.inv %*% A))) {
    list(basis=matrix(0, nrow=ncol(A), ncol=0),
         translate=as.vector(A.inv %*% b))
  } else {
    the.qr <- qr(diag(n) - A.inv %*% A)
    list(basis=qr.Q(the.qr)[, 1:the.qr$rank, drop=FALSE],
         translate=as.vector(A.inv %*% b))
  }
}



# Generate a projection matrix that transforms an (n-1) dimensional vector in
# homogeneous coordinate representation to an n-dimensional weight vector.
createTransform <- function(basis, inverse=FALSE, keepHomogeneous=inverse) {
  # add one extra element to vectors in each basis (homogeneous coordinate
  # representation)
  translate <- if (inverse == FALSE) basis$translate else -basis$translate
  n <- length(translate)
  basis <- basis$basis
  basis <- rbind(cbind(basis, rep(0, n)), c(rep(0, ncol(basis)), 1))
  
  # create translation matrix (using homogenous coordinates)
  translation <- rbind(cbind(diag(n), translate), c(rep(0, n), 1))
  
  # homogeneous coordinate elimination
  nh <- if (inverse == FALSE) { nrow(basis) } else { ncol(basis) }
  elimHom <- if (keepHomogeneous) {
    diag(nh)
  } else {
    cbind(diag(nh - 1), rep(0, nh - 1))
  }
  
  transform <- if (inverse == FALSE) {
    # successively apply basis transformation and translation
    elimHom %*% translation %*% basis
  } else {
    # successively apply translation and basis transformation
    elimHom %*% t(basis) %*% translation
  }
  
  unname(transform)
}


transformConstraints <- function(transform, constr) {
  list(
    constr = constr$constr %*% transform,
    dir = constr$dir,
    rhs = constr$rhs)
}



checkPolytope <- function(x0, constr, homogeneous, transform) {
  n <- if (!is.null(x0)) length(x0) else ncol(constr$constr)
  m <- nrow(constr$constr)
  
  # Verify preconditions
  stopifnot(n > homogeneous)
  stopifnot(n == ncol(constr$constr))
  stopifnot(m == length(constr$rhs))
  stopifnot(constr$dir == "<=")
  
  if (homogeneous) { # Change to homogeneous coordinates
    stopifnot(x0[n] == 1.0)
    list(n = n - 1,
         m = m,
         x0 = if (is.null(x0)) x0 else x0[1:(n - 1)],
         constr = list(constr = constr$constr[ , 1:(n - 1), drop=FALSE],
                       rhs = constr$rhs - constr$constr[ , n, drop=TRUE],
                       dir = constr$dir),
         transform = function(samples) {
           if (!is.null(transform)) {
             mat <- samples %*% t(transform[ , 1:(n - 1), drop=FALSE])
             t(t(mat) + transform[ , n, drop=TRUE])
           } else {
             cbind(samples, 1)
           }
         },
         xN = function(samples) {
           c(samples[nrow(samples), , drop=TRUE], 1)
         })
  } else {
    list(n = n,
         m = m,
         x0 = x0,
         constr = constr,
         transform = function(samples) {
           if (!is.null(transform)) {
             samples %*% t(transform)
           } else {
             samples
           }
         },
         xN = function(samples) {
           samples[nrow(samples), , drop=TRUE]
         })
  }
}


har <- function(x0, constr, N, thin=1, homogeneous=FALSE, transform=NULL) {
  stopifnot(N %% thin == 0)
  args <- checkPolytope(x0, constr, homogeneous, transform)
  
  rval <- .Call(hitandrun_har, args$x0, args$constr$constr, args$constr$rhs, N, thin)
  
  list(samples=args$transform(rval),
       xN=args$xN(rval))
}
