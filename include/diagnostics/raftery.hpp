// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RAFTERY_HPP
#define RAFTERY_HPP

template <typename NT>
NT fix(NT x)
{
    if (x > 0.0) return std::floor(x);

    return std::ceil(x);
}

#include "raftery_subroutines/empquant.hpp"
#include "raftery_subroutines/indtest.hpp"
#include "raftery_subroutines/mctest.hpp"
#include "raftery_subroutines/mcest.hpp"
#include "raftery_subroutines/thin.hpp"
#include "raftery_subroutines/ppnd.hpp"



template <typename VT, typename MT, typename NT>
std::pair<MT,VT> perform_raftery(MT samples, NT q, NT r, NT s)
{
    typedef Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> MTint;
    typedef Eigen::Matrix<int,Eigen::Dynamic,1> VTint;

    //MT samples = samples2.transpose();
    unsigned int n = samples.rows(), d = samples.cols(), kthin, kmind;
    //std::cout<<"n = "<<n<<", d = "<<d<<std::endl;
    MT results(d,5);
    VT res(d);
    MTint work = MTint::Zero(n,d); 
    VTint tmp = VTint::Zero(n);
    std::pair<int,VTint> xy;
    std::pair<NT,NT> g2bic;

    NT cutpt, alpha, beta, g2, bic, epss;
    int tcnt;

    for (int i =0; i<d; i++)
    {
        cutpt = empquant<VT>(samples, q);
        for (int j = 0; j<n; j++)
        {
            for (int k = 0; k<d; k++)
            {
                if (samples(j,k) <= cutpt) work(j,k) = 1;
            }
        }
        kthin = 1; bic = 1.0; epss = 0.001;

        while(bic > 0.0) 
        {
            xy = thin<VTint>(work, n, kthin);
            tcnt = xy.first;
            tmp = xy.second;
            g2bic = mctest<MTint, NT>(tmp, tcnt);
            g2 = g2bic.first;
            bic = g2bic.second;
            kthin++;
            if (kthin > n/2) {
                break;
            }
        }
        
        kthin--;
        g2bic = mcest<MTint, NT>(tmp, tcnt);
        alpha = g2bic.first;
        beta = g2bic.second;
        kmind = kthin;
        g2bic = indtest<MTint, NT>(tmp, tcnt);
        g2 = g2bic.first;
        bic = g2bic.second;

        while (bic> 0.0)
        {
            xy = thin<VTint>(work, n, kmind);
            g2bic = indtest<MTint, NT>(tmp, tcnt);
            g2 = g2bic.first;
            bic = g2bic.second;
            kmind++;
            if (kmind > n) {
                break;
            }
        }

        NT psum = alpha + beta;
        NT tmp1  = std::log(psum*epss/std::max(alpha,beta)) / std::log(std::abs(1.0 - psum));
        NT nburn = fix((tmp1+1.0)*NT(kthin));
        NT phi   = ppnd((s+1.0)/2.0);
        NT tmp2  = (2.0 - psum)*alpha*beta*(phi*phi)/(psum*psum*psum * (r*r));
        NT nprec = fix(tmp2+1.0)*kthin; 
        NT nmin  = fix(((1.0-q)*q*(phi*phi)/(r*r))+1.0);
        NT irl   = (nburn + nprec) / nmin; 
        NT kind  = std::max(fix(irl+1.0), NT(kmind));

        results(i,0) = NT(kthin);
        results(i,1) = NT(nburn);
        results(i,2) = kind;
        results(i,3) = NT(nburn)+nprec;
        results(i,4) = nmin;

        res(i) = irl;


    }

    return std::pair<MT,VT>(results, res);

}


#endif
