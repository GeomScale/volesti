**VolEsti** is a `C++` library for volume approximation and sampling of convex bodies (*e.g.* polytopes) with an `R` and limited `python` interface. **VolEsti** is part of the [GeomScale](https://geomscale.github.io) project.

[![CRAN status](https://www.r-pkg.org/badges/version/volesti)](https://cran.r-project.org/package=volesti)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/volesti)](https://cran.r-project.org/package=volesti)
![CRAN/METACRAN](https://img.shields.io/cran/l/volesti)
[![Chat](https://badges.gitter.im/boostorg/geometry.png)](https://gitter.im/GeomScale/community?utm_source=share-link&utm_medium=link&utm_campaign=share-link)

### Test results

[![CRAN check](https://cranchecks.info/badges/worst/volesti)](https://cran.r-project.org/web/checks/check_results_volesti.html)
[![CircleCI master](https://circleci.com/gh/GeomScale/volume_approximation/tree/master.svg?style=shield)](https://circleci.com/gh/GeomScale/volume_approximation/tree/master)
[![CircleCI master](https://circleci.com/gh/GeomScale/volume_approximation/tree/develop.svg?style=shield)](https://circleci.com/gh/GeomScale/volume_approximation/tree/develop)

###  Documentation

* [Using the R interface](doc/r_interface.md)
* [Using the C++ Interface](doc/cpp_interface.md)
* [Wikipage with tutorials and demos](https://github.com/GeomScale/volume_approximation/wiki)
* [Tutorial given to PyData meetup](https://vissarion.github.io/tutorials/volesti_tutorial_pydata.html)

### Credits

Copyright (c) 2012-2020 Vissarion Fisikopoulos  
Copyright (c) 2018-2020 Apostolos Chalkis  
Copyright (c) 2020-     Marios Papachristou

You may redistribute or modify the software under the GNU Lesser General Public License as published by Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.  

Main development started in 2012 by Vissarion Fisikopoulos while he was affiliated with University of Athens (UoA, Greece), and then with University of Brussels (ULB, Belgium). The main sampling and volume algorithms were developed at that time. 

Later, Chalkis Apostolos affiliated with University of Athens (UoA, Greece) and mentored by V.F. and partially supported by GSoC'2018, 2019 programs enhanced volesti with more algorithms for sampling and volume as well as an `R` interface.       

We acknowledge several contributions by the open-source community, most notably a `python` interface by Pedro Zuidberg Dos Martires affiliated with KU Leuven.

### Publications

1. I.Z. Emiris and V. Fisikopoulos, *Efficient random-walk methods for approximating polytope volume*, In Proc. ACM Symposium on Computational Geometry, Kyoto, Japan, p.318-325, 2014.  
2. I.Z. Emiris and V. Fisikopoulos, *Practical polytope volume approximation*, ACM Transactions on Mathematical Software, vol 44, issue 4, 2018.
3. L. Cales, A. Chalkis, I.Z. Emiris, V. Fisikopoulos, *Practical volume computation of structured convex bodies, and an application to modeling portfolio dependencies and financial crises*, Proc. of Symposium on Computational Geometry, Budapest, Hungary, 2018.  

