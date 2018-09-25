// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#ifndef INITIALIZATION_H
#define INITIALIZATION_H

typedef double NT;
typedef Cartesian<NT>    Kernel;
typedef typename Kernel::Point    Point;
typedef boost::mt19937    RNGType;
typedef HPolytope<Point> Hpolytope;
typedef VPolytope<Point, RNGType > Vpolytope;
typedef Zonotope<Point> Zonotope;
typedef copula_ellipsoid<Point> CopEll;

#endif
