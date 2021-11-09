
#' @export
is_child <- function(S_indices, node_indices) {
  
  ischild = FALSE
  
  n = length(S_indices)
  
  for (i in 1:n) {
    indx = S_indices[i]
  
    if (sum(indx %in% node_indices) > 0) {
      ischild = TRUE
      return(ischild)
    }
  }
  return(ischild)
}

#' @export
is_leaf <- function(node, X, S, S_Vindices, V, A, b, x0, lb, ub, nu) {
  
  res = list()
  isleaf = TRUE
  S_ind_to_esti = c()
  number_of_leaves = 0
  ratios = c()
  all_ratios = c()
  S_Vindices_leaves = matrix(list(), 0, 1)
  counter = 1
  n = length(S)
  
  for (i in 1:n) {
    Sind = S_Vindices[[i]]
    conv = FALSE
    if (is_child(S_Vindices[[i]], node$V_indices)) {
      check_res = check_convergence_interface(A, b, V, Sind, x0, X, nu, lb, ub, TRUE)
      all_ratios = c(all_ratios, check_res$ratio)
    }
    if (!check_res$conv) {
      isleaf = FALSE
    } else {
      number_of_leaves = number_of_leaves + 1
      S_ind_to_esti = c(S_ind_to_esti, i)
      ratios = c(ratios, check_res$ratio)
      S_Vindices_leaves[[counter]] = S_Vindices[[i]]
      counter = counter + 1
    }
  }
  
  res$is_leaf = isleaf
  res$S_ind_to_esti = S_ind_to_esti
  res$S_Vindices_leaves = S_Vindices_leaves
  res$number_of_leaves = number_of_leaves
  res$ratios = ratios
  res$all_ratios = all_ratios
  
  return(res)
}

#' @export
get_next_child <- function(node, X, cmin, cmax, S_Vindices, A, b, V, x0, lb, ub, nu) {
  
  cmin_temp = cmin
  cmax_temp = cmax
  converged = FALSE
  iter = 1
  last_round = FALSE
  
  while(!converged & !last_round) {
  
    if (iter == 15) {
      last_round = TRUE
    }
    
    cmed = (cmax_temp + cmin_temp) / 2
    Vmed = cmed * V
    b_med = cmed * b
    
    find_components_res = find_components_new_vertices(Vmed, x0)
    S_med = find_components_res$S
    S_Vindices_med = find_components_res$S_Vindices
    
    important_res = keep_important_components(S_Vindices_med, S_Vindices, Vmed)
    S_med = important_res$S
    S_Vindices_med = important_res$S_Vindices
    
    num_of_components = length(S_med)
    ratios = c()
    converges = c()
    Xs = matrix(, nrow=length(x0), ncol=0)
    too_small_ratios = c()
    
    for (i in 1:num_of_components) {
      
      Sind = S_Vindices_med[[i]]
      check_res = check_convergence_interface(A, b_med, Vmed, Sind, x0, X, nu, lb, ub, last_round)

      too_small_ratios = c(too_small_ratios, check_res$too_few)
      converges = c(converges, check_res$conv)
      ratios = c(ratios, check_res$ratio)
      Xs = cbind(Xs, check_res$x)
    }
    
    indx = which.min(ratios)
    
    if (too_small_ratios[indx]) {
      cmin_temp = cmed
      iter = iter+1
      next
    }
    if (converges[indx]) {
      converged = TRUE
      break
    } else {
      cmax_temp = cmed
    }
    iter = iter+1
  }
  
  if ((converged && (ratios[indx] < 0.9)) || (!converged && last_round && (ratios[indx] > 0.03))) {
    next_node = GenerateNode(cmed, Xs[,indx], S_med[[indx]], S_Vindices_med[[indx]], ratios[indx])
  } else {
    next_node = list()
  }
  
  return(next_node)
}


#' @export
build_tree <- function(node, X, cmin, cmax, A, b, S, S_Vindices, V, x0, lb, ub, nu, N, W) {
  
  tree = list()
  if(!is_father_of_a_leaf(node, S_Vindices)) {
    tree = list()
    return(tree)
  }
  
  while (TRUE) {
    
    if (length(S) == 0) {
      tree = node
      return(tree)
    }
    
    next_node = get_next_child(node, X, cmin, cmax, S_Vindices, A, b, V, x0, lb, ub, nu)
    
    if (length(next_node) == 0) {
      is_leaf_res = is_leaf(node, X, S, S_Vindices, V, A, b, x0, lb, ub, nu)
    
      indx = which.min(is_leaf_res$all_ratios)
      remove_leaf_res = remove_leaf(S, S_Vindices, indx)
      S = remove_leaf_res$S
      S_Vindices = remove_leaf_res$S_Vindices
      next
    }
    
    Y = sample_component(A, next_node$c*b, next_node$x0, N, W, x0)
    
    is_leaf_res = is_leaf(next_node, Y, S, S_Vindices, V, A, b, x0, lb, ub, nu)
    
    if (is_leaf_res$number_of_leaves > 0) {
      next_node$S_ind_to_esti = is_leaf_res$S_ind_to_esti
      next_node$number_of_leaves = is_leaf_res$number_of_leaves
      next_node$ratio_of_leaves = is_leaf_res$ratios
      next_node$S_Vindices_leaves = is_leaf_res$S_Vindices_leaves
    
      if (is_leaf_res$is_leaf) {
        next_node$IsLeaf = TRUE
      }
      
      remove_leaf_res = remove_leaf(S, S_Vindices, is_leaf_res$S_ind_to_esti)
      
      S = remove_leaf_res$S
      S_Vindices = remove_leaf_res$S_Vindices

      if (length(S) == 0) {
        node$children[[length(node$children)+1]] = next_node
        tree = node
        return(tree)
      } else {
        node$children[[length(node$children)+1]] = next_node
      }
    } else {
      node$children[[length(node$children)+1]] = next_node
    }
    
    cmax = next_node$c
    if (is_leaf_res$is_leaf) {
      next
    } else {
      if(is_father_of_a_leaf(next_node, S_Vindices)) {
        
        leaves_from_node_res = get_leaves_from_node(S, S_Vindices, next_node)
        S2 = leaves_from_node_res$S
        S_Vindices2 = leaves_from_node_res$S_Vindices
        
        subtree = build_tree(next_node, Y, cmin, cmax, A, b, S2, S_Vindices2, V, x0, lb, ub, nu, N, W)

        node$children[[length(node$children)]] = subtree
        
        remove_leaf_res = remove_leaves_from_node(S, S_Vindices, next_node)
        
        S = remove_leaf_res$S
        S_Vindices = remove_leaf_res$S_Vindices
      }
    }
  }
  
  return(tree)
}



#' @export
get_tree_of_components <- function(V, A, b, x0, lb, ub, nu, N, W) {
  
  res = list()
  single_node = FALSE
  empty = FALSE
  cmax = get_c_upper_bound(A, b, x0)
  cmin = 1
  components_res = find_components_new_vertices(V, x0)
  S = components_res[[1]]
  S_Vindices = components_res[[2]]
  if(length(components_res) > 2){
    V_ind_out = components_res[[3]]
  } else {
    V_ind_out = c()
  }
  
  if (length(S_Vindices) == 0) {
    empty = TRUE
    res$single_node = single_node
    tree = list()
    res$empty = empty
    res$tree = tree
    return(res)
  }
  #print(length(S))
  #print(length(S_Vindices[[1]]))
  if (length(S) == 1){
    single_node = TRUE
    interior_res = compute_interior_point_in_node_eff(S_Vindices[[1]], V_ind_out, V, x0, A, b)
    if (interior_res$found) {
      tree = GenerateNode(1, interior_res$x, S[[1]], S_Vindices[[1]], NaN)
      res$single_node = single_node
      res$tree = tree
      res$empty = FALSE
      return(res)
    } else {
      #X = boundary_randsphere(d, N) + kronecker(matrix(1, 1, N), matrix(x0, ncol = 1))
      fast_res = get_fast_interior_point(A, b, x0, V, S_Vindices[[1]], cmin, cmax)
      tree = GenerateNode(1, fast_res, S[[1]], S_Vindices[[1]], NaN)
      res$single_node = single_node
      res$tree = tree
      res$empty = FALSE
      return(res)
    }
    #N=1200
    #ub=0.7
    #print(N)
    #print(ub)
  }
  
  #N=1200
  #print(N)
  
  n = dim(V)[2]
  d = length(x0)
  
  X = boundary_randsphere(d, N) + kronecker(matrix(1, 1, N), matrix(x0, ncol = 1))
  
  node = GenerateNode(cmax, X[,1], cmax*V, 1:n, 1)
  
  is_leaf_res = is_leaf(node, X, S, S_Vindices, V, A, b, x0, lb, ub, nu)
  node$isLeaf = is_leaf_res$is_leaf
  
  if (is_leaf_res$number_of_leaves > 0) {
    node$S_ind_to_esti = is_leaf_res$S_ind_to_esti
    node$number_of_leaves = is_leaf_res$number_of_leaves
    node$ratio_of_leaves = is_leaf_res$ratios
    node$S_Vindices_leaves = is_leaf_res$S_Vindices_leaves
  
    remove_leaf_res = remove_leaf(S, S_Vindices, S_ind_to_esti)
    S = remove_leaf_res$S
    S_Vindices = remove_leaf_res$S_Vindices
  }
  
  if (!node$isLeaf) {
    tree = build_tree(node, X, cmin, cmax, A, b, S, S_Vindices, V, x0, lb, ub, nu, N, W)
    num_leaves = count_num_of_leaves(tree, 0)
    
    if (num_leaves == 1) {
      single_node = TRUE
    }
  } else {
    tree = node
  }
  
  res$single_node = single_node
  res$tree = tree
  res$empty = FALSE
  return(res)
}


#' @export
get_samples_from_leaf <- function(node, A, b, x0, V, N, W, target_psrf) {
  
  while(TRUE) {
    if (node$number_of_leaves > 0) {
      break
    }
    node = node$children[[1]]
  }
  
  X = sample_component(A, node$c*b, node$x0, 1000, 10, x0)
  
  cols = node$S_Vindices_leaves[[1]]
  
  x = return_first_inside_interface(A, b, V, cols, x0, X)
  
  node = GenerateNode(1, x, V, cols, node$ratio_of_leaves[[1]])
  
  #samples = sample_component(A, b, node$x0, N, W, x0)
  samples = sample_with_psrf(node, A, b, x0, N, W, psrf_target)
  
  return(samples)
}


#' @export
get_samples_from_tree <- function(tree, single_node, A, b, x0, V, N, W, psrf_target) {
  
  samples = matrix(,0,0)
  
  if (single_node) {
    if (length(tree$children) == 0) {
      node = tree
      #samples = sample_component(A, b, node$x0, N, W, x0)
      samples = sample_with_psrf(node, A, b, x0, N, W, psrf_target)
    } else {
      samples = get_samples_from_leaf(tree, A, b, x0, V, N, W, psrf_target)
      #samples = sample_with_psrf(tree, A, b, x0, N, W, psrf_target)
    }
  } else {
    print(paste0('No single_node'))
  }
  return(samples)
}


#' @export
sample_with_psrf <- function(node, A, b, x0, N, W, psrf_target) {
  
  samples = sample_component_psrf(A,b,node$x0,2*N,W,x0,psrf_target)
  #print(psrf_univariate(samples))
  
  if (dim(samples)[2] > N) {
    indx <- sample(1:dim(samples)[2], N)
    samples = samples[,indx]
  }
  
  return(samples)
}


#' @export
sample_with_psrf_2 <- function(node, A, b, x0, N, W, psrf_target) {
  
  #samples = sample_component_psrf(A,b,node$x0,2*N,20,x0,psrf_target)
  #print(psrf_univariate(samples))
  MAX_ITER = 1000
  iter = 1
  psrf_val = 2 * psrf_target
  d = dim(A)[2]
  samples = matrix(,d,0)
  
  while(iter <= MAX_ITER && psrf_val > psrf_target) {
    
    X = sample_component(A, b, node$x0, 2*N, W, x0)
    samples = cbind(samples, X)
    psrf_val = max(psrf_univariate(samples))
    print(psrf_val)
    iter = iter + 1
  }
  
  if (dim(samples)[2] > N) {
    indx <- sample(1:dim(samples)[2], N)
    samples = samples[,indx]
  }
  
  return(samples)
}




