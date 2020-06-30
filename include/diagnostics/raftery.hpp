// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RAFTERY_HPP
#define RAFTERY_HPP

template <typename VT, typename MT, typename NT>
std::pair<MT,VT> perform_raftery(MT samples2, NT q, NT r, NT s)
{
    MT samples = samples2.transpose();
    unsigned int n = samples.rows(), d = samples.cols(), kthin, kmind;
    MT work(n,d), results(d,4);
    VT tmp = VT::Zero(n), bic, epss, res(d);
    std::pair<NT,VT> xy;
    std::pair<NT,NT> g2bic;

    NT cutpt, alpha, beta, tcnt, g2;

    for (int i =0; i<d; i++)
    {
        cutpt = empquant(samples, q);
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
            xy = thin(work, n, kthin);
            tcnt = xy.first;
            tmp = xy.second;

            g2bic = mctest(tmp, tcnt);
            g2 = g2biv.first;
            bic = g2bic.second;
            kthin++;
        }
        
        kthin--:
        g2bic = mcest(tmp, tcnt);
        alpha = g2bic.first;
        beta = g2bic.second;
        kmind = kthin;
        g2bic = indtest(tmp, tcnt);
        g2 = g2bic.first;
        bic = g2bic.second;

        while (bic> 0.0)
        {
            xy = thin(woek, n, kmind);
            g2bic = indtest(tmp, tcnt);
            g2 = g2bic.first;
            bic = g2bic.second;
            kmind++;
        }

        NT psum = alpha + beta;
        NT tmp1  = std::log(psum*epss/std::max(alpha,beta)) / std::log(std::abs(1.0 - psum));
        NT nburn = std::floor((tmp1+1.0)*NT(kthin));
        NT phi   = ppnd((s+1.0)/2.0);
        NT tmp2  = (2.0 - psum)*alpha*beta*(phi*phi)/(psum*psum*psum * (r*r));
        NT nprec = std::floor(tmp2+1.0)*kthin; 
        NT nmin  = std::floor(((1.0-q)*q*(phi*phi)/(r*r))+1.0);
        NT irl   = (nburn + nprec) / nmin; 
        NT kind  = std::max(std::floor(irl+1.0), NT(kmind));

        results(i,0) = NT(kthin);
        results(i,1) = NT(nburn);
        results(i,2) = NT(nburn)+nprec);
        results(i,3) = nmin;

        res(i) = irl;


    }

    return std::pair<MT,VT>(results, res);

}


#endif
