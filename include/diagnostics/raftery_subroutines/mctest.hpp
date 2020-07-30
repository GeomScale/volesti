// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

//Based on Matlab version of coda package in http://www.spatial-econometrics.com/gibbs/

#ifndef MCTEST_HPP
#define MCTEST_HPP


template<typename MT, typename NT, typename VT>
std::pair<NT,NT> mctest(VT const& d, int const& n)
{
   MT m1 = MT::Zero(2,2), m2 = MT::Zero(2,2);
   NT g2 = 0.0, bic = 0.0, fitted;
   int i1, i2, i3, t1, t2, t3, t4, focus;

   for (int i=2; i<n; i++)
   {
      i1 = d(i-2)+1; i2 = d(i-1)+1; i3 = d(i)+1;
      if (i3 == 1) m1(i1-1, i2-1) = m1(i1-1, i2-1) + 1;
      if (i3 == 2) m2(i1-1, i2-1) = m2(i1-1, i2-1) + 1;
   }
   for (i1 = 0; i1 < 2; i1++) {
      for (i2 = 0; i2 < 2; i2++) {
         for (i3 = 0; i3 < 2; i3++) {
            if (i3 == 0){
               if (m1(i1, i2) != 0){
                  t1 = m1(i1, i2)+m2(i1, i2); 
                  t2 = m1(0, i2)+m1(1, i2);
                  t3 = m1(0, i2)+m2(0, i2);
                  t4 = m1(1, i2)+m2(1, i2);
                  fitted = (NT(t1) * NT(t2)) / (NT(t3) + NT(t4)); 
                  focus = m1(i1,i2); 
                  g2 = g2 + std::log(NT(focus) / fitted) * NT(focus); 
               }
            }
            if (i3 == 1) {
               if (m2(i1, i2) != 0) {
                  t1 = m1(i1, i2) + m2(i1, i2); 
                  t2 = m2(0, i2) + m2(1, i2);
                  t3 = m1(0, i2) + m2(0, i2);
                  t4 = m1(1, i2) + m2(1, i2);
                  fitted = (NT(t1) * NT(t2)) / (NT(t3) + NT(t4)); 
                  focus = m2(i1, i2); 
                  g2 = g2 + std::log(NT(focus) / fitted) * NT(focus);
               }
            }
         }
      }
   }
   g2 = g2 * 2.0;
   bic = g2 - std::log(NT(n) - 2.0) * 2.0;
   return std::pair<NT,NT>(g2, bic);
}


#endif
