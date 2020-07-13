// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef GET_FULL_DIMENSIONAL_POLYTOPE
#define GET_FULL_DIMENSIONAL_POLYTOPE


template <typename H_polytope, typename MT, typename VT>
std::pair<H_polytope, std::pair<MT, VT> > get_full_dimensional_polytope(MT A, VT b, MT Aeq, VT beq, bool slow = false)
{
    typedef typename H_polytope::NT NT;

    VT p = Aeq.colPivHouseholderQr().solve(beq);
    int r = Aeq.rows(), d = A.cols();
    MT V;
    VT s;
    
    //std::cout<<"hello"<<std::endl;
    
    if (slow){
        Eigen::JacobiSVD<MT> svd(Aeq, Eigen::ComputeFullV);
        V = svd.matrixV();
        s = svd.singularValues();
    } else {
        //std::cout<<"hello2"<<std::endl;
        Eigen::BDCSVD<MT> bdcSvd(Aeq, Eigen::ComputeFullV);
        V = bdcSvd.matrixV();
        s = bdcSvd.singularValues();
        //std::cout<<"hello3"<<std::endl;
    }

    const NT e = 0.0000000000001;
    int r_count = 0;
    for (int i=0 ; i<s.rows() ; i++) {
        if (std::abs(s(i)) <= e){
            r_count++;
            //N.conservativeResize(N.rows(), N.cols()+1);
            //N.col(N.cols()-1) = V.col(i);
        }
    }
    //std::cout<<"hello4"<<std::endl;
    MT N(d, d - r + r_count);
    N = V.block(0, r - r_count, d, d - r + r_count);
    //std::cout<<"hello5"<<std::endl;
    
    //for (int i = s.rows(); i<V.cols(); i++){
        //N.conservativeResize(N.rows(), N.cols()+1);
        //N.col(N.cols()-1) = V.col(i);
    //}

    //std::cout<<"N = "<<N.rows()<<"N.cols() = "<<N.cols()<<std::endl;
    b = b - A * p;
    //MT A2 = A * N;

    H_polytope HP;
    HP.init(N.cols(), A * N, b);

    return std::pair<H_polytope, std::pair<MT, VT> >(HP, std::pair<MT,VT>(N, p));
}

#endif
