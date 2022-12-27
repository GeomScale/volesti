sample_from_last_polytope <- function(mmcs_output_list, N) {

    inner_ball = get_max_inner_ball(mmcs_output_list$HP_rounded$A, mmcs_output_list$HP_rounded$b)
    more_samples = sample_points(mmcs_output_list$HP_rounded, n = N, 
                                 random_walk = list("walk" = "aBiW", "walk_length" = 1,
                                                    "starting_point"=inner_ball$center, 
                                                    "L" = 6*sqrt(mmcs_output_list$HP_rounded$dimension)*inner_ball$radius))

    ## map the points back to P0
    samples_in_P0 = mmcs_output_list$T %*% more_samples + 
        kronecker(matrix(1, 1, N), matrix(mmcs_output_list$T_shift, ncol = 1))

    ## compute the portfolios
    more_random_portfolios = mmcs_output_list$N %*% samples_in_P0 + 
        kronecker(matrix(1, 1, N), matrix(mmcs_output_list$N_shift, ncol = 1))

    return(more_random_portfolios)
}