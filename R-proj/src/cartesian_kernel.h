#include <Rcpp.h>
#include "point.h"


#ifndef CARTESIAN_KERNEL_H
#define CARTESIAN_KERNEL_H



template <typename K>
class Cartesian
{
public:
  typedef Cartesian<K> Self;
  typedef K                    RT;
  typedef point<Self>              Point;
    
};

#endif

