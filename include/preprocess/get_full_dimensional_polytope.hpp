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

   //VT p = Aeq.colPivHouseholderQr().solve(beq), s;
   MT AA = Aeq;
   VT p = AA.colPivHouseholderQr().solve(beq), s;
   
   int r = Aeq.rows(), d = A.cols();
   MT V;
    
   if (slow){ //typically better when dimension <= 15
      Eigen::JacobiSVD<MT> svd(Aeq, Eigen::ComputeFullV);
      V = svd.matrixV();
      s = svd.singularValues();
   } else {
      Eigen::BDCSVD<MT> bdcSvd(Aeq, Eigen::ComputeFullV);
      V = bdcSvd.matrixV();
      s = bdcSvd.singularValues();
   }

   const NT e = 0.0000000000001;
   int r_count = 0;
   
   if (d > r) {
      r_count = d - r;
   }

   for (int i=0 ; i<s.rows() ; i++) { // count zero singular values
      if (std::abs(s(i)) <= e){
         r_count++;
      }
   }

   //MT N(d, d - r + r_count); // the null space is the columns of V that correspond to zero singular values
   //N = V.block(0, r - r_count, d, d - r + r_count);
   
   std::cout<<"r_count = "<<r_count<<"\n"<<std::endl;


   MT N(d, r_count); // the null space is the columns of V that correspond to zero singular values
   N = V.block(0, d - r_count, d, r_count);

   b = b - A * p;

   H_polytope HP(N.cols(), A * N, b);

   return std::pair<H_polytope, std::pair<MT, VT> >(HP, std::pair<MT,VT>(N, p));
}

#endif
