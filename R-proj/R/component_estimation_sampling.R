
#' export
get_target_estimation_errors <- function(node, error, depth = 1) {
  
  has_leaf = 0
  if (node$number_of_leaves > 0) {
    node$error_to_leaf = error / (sqrt(depth+1-1))
    has_leaf = 1
  }
  
  if(length(node$children) == 0) {
    node$error_to_estimate = error / (sqrt(depth + has_leaf -1))
    node$max_depth = depth + has_leaf
  } else {
    n = length(node$children)
    for (i in 1:n) {
      node2 = node$children[[i]]
      child = get_target_estimation_errors(node2, error, depth+1)
      if (child$max_depth > node$max_depth) {
        node$max_depth = child$max_depth
        node$error_to_estimate = error / sqrt(node$max_depth - 1)
      }
      node$children[[i]] = child
    }
  }
  
  node$depth = depth
  return(node)
}


#' export
compute_node_weights <- function(node, A, b, x0, V, win_len, N, ratio_estimated = 1, all_leaves_info = matrix(list(), 0, 1)) {
  
  d = length(x0)
  leaves = list()
  leaves$ratio_leaves = c()
  leaves$Xs = matrix(,d,0)
  leaves$leaf_verts = matrix(list(), 0, 1)
  storing = FALSE
  
  if (node$number_of_leaves > 0) {
    for (i in 1:node$number_of_leaves) {
      esti_res = estimate_component(A, b, node$x0, 1, win_len, x0, V, node$c, 1, node$error_to_leaf, node$ratio_of_leaves[i], N, storing)
      leaves$ratio_leaves = c(leaves$ratio_leaves, esti_res$ratio * ratio_estimated)
      leaves$Xs = cbind(leaves$Xs, esti_res$x)
      leaves$leaf_verts[[i]] = node$S_Vindices_leaves[[i]]
    }
    #leaves$ratios = ratio_leaves
    #leaves$leaf_verts = leaf_verts
    all_leaves_info[[length(all_leaves_info) + 1]] = leaves
  }
  
  if(length(node$children) > 0) {
    n = length(node$children)
    for (i in 1:n) {
      child = node$children[[i]]
      esti_res = estimate_component(A, b, node$x0, 1, win_len, x0, V, node$c, child$c, child$error_to_estimate, child$ratio, N, storing)
      all_leaves_info = compute_node_weights(child, A, b, x0, V, win_len, N, esti_res$ratio, all_leaves_info)
    }
  }
  
  return(all_leaves_info)
}


#' export
sample_from_leaves <- function(all_leaves_info, A, b, x0, N, W, psrf_target) {
  
  n = length(all_leaves_info)
  print(n)
  d = length(x0)
  samples = matrix(list(), 0, 1)
  final_samples = matrix(0, d, N)
  sample_lens = c()
  ratios = c()
  
  for (i in 1:n) {
    leaf = all_leaves_info[[i]]
    for (j in 1:length(leaf$ratio_leaves)) {
      ratios = c(ratios, leaf$ratio_leaves[j])
      x = leaf$Xs[,j]
      X = sample_component_psrf(A,b,x,2*N,W,x0,psrf_target)
      sample_lens = c(sample_lens, dim(X)[2])
      samples[[length(samples) + 1]] = X
    }
  }
  
  rel_ratios = ratios / sum(ratios)
  Us = runif(N, 0, 1)
  cum_sum_ratios = cumsum(rel_ratios)
  
  for (i in 1:N) {
    indx = which(Us[i] < cum_sum_ratios)[1]
    indx_col <- sample(1:sample_lens[indx], 1)
    final_samples[, i] = (samples[[indx]])[, indx_col]
  }
  
  return(final_samples)
}



