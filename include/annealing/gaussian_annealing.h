#ifndef GAUSSIAN_ANNEALING_H
#define GAUSSIAN_ANNEALING_H

#include <complex>

template <class T1>
int get_first_gaussian(T1 &K, NT &a0, NT frac, vars var){

    int m=K.num_of_hyperplanes(), dim=var.n;
    NT sum, lower=0.0, upper=1.0;
    std::vector dists(m,0);
    for(int i=0; i<m; i++){
        sum=0.0;
        for(int j=1; j<dim+1; j++){
            sum+=P.get_coeff(i,j)*P.get_coeff(i,j);
        }
        dists[i]=P.get_coeff(i,0)/std::sqrt(sum);
    }
    return 1;
}


#endif
