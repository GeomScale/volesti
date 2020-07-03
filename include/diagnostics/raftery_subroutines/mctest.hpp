// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef MCTEST_HPP
#define MCTEST_HPP


template<typename NT, typename VT>
std::pair<NT,NT> mctest(VT d, int n)
{
     m1 = VT::Zero(2,2); m2 = VT::Zero(2,2);
     NT g2 = 0.0, bic = 0.0;
     int i1, i2, i3, t1, t2, t3, t4, fitted, focus;

 for (int i=2; i<n; i++)//=3:n;      % count states
 {
    i1 = d(i-2,0)+1; i2 = d(i-1,0)+1; i3 = d(i,0)+1;
    if (i3 == 1) m1(i1,i2) = m1(i1,i2) + 1;// end;
    if (i3 == 2) m2(i1,i2) = m2(i1,i2) + 1;// end;
 }//       % end of for i
 for (i1 = 0; i1<2; i1++){// 1:2;
  for (i2 = 0; i2<2; i2++){//1:2;
   for (i3 = 0; i3<2; i3++){//} 1:2;
   if (i3 == 0){
    if (m1(i1,i2) != 0){
    t1 = m1(i1,i2)+m2(i1,i2); t2 = m1(1,i2)+m1(2,i2);
    t3 = m1(1,i2)+m2(1,i2); t4 = m1(2,i2)+m2(2,i2);
    fitted = (t1*t2)/(t3+t4); focus = m1(i1,i2); 
    g2 = g2 + std::log(NT(focus)/NT(fitted)) * NT(focus); 
    }//end;      % end of if m1   
   }//end;       % end of if i3
   if (i3 == 1) {
    if (m2(i1,i2) != 0) {
    t1 = m1(i1,i2) + m2(i1,i2); t2 = m2(1,i2) + m2(2,i2);
    t3 = m1(1,i2) + m2(1,i2); t4 = m1(2,i2) + m2(2,i2);
    fitted = (t1*t2)/(t3+t4); focus = m2(i1,i2); 
    g2 = g2 + std::log(NT(focus) / NT(fitted)) * NT(focus);
    }//end;      % end of if m2  
   }//end;        % end of if i3
   }//end;       % end of for i3 
  }//end;        % end of for i2 
 }//end;         % end of for i1
 g2 = g2*2.0; bic = g2 - std::log(NT(n) - 2.0)*2.0;

 return std::pair<NT,NT>(g2, bic);
}



#endif
