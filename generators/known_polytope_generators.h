// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef KNOWN_POLYTOPE_GENERATORS_H
#define KNOWN_POLYTOPE_GENERATORS_H

#include <exception>

template <typename Polytope>
Polytope gen_cube(const unsigned int &dim, const bool &Vpoly) {

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;
    MT A;
    VT b;
    unsigned int m;

    if (!Vpoly) {

        A.resize(2 * dim, dim);
        b.resize(2 * dim);
        for (unsigned int i = 0; i < dim; ++i) {
            b(i) = 1.0;
            for (unsigned int j = 0; j < dim; ++j) {
                if (i == j) {
                    A(i, j) = 1.0;
                } else {
                    A(i, j) = 0.0;
                }
            }
        }
        for (unsigned int i = 0; i < dim; ++i) {
            b(i + dim) = 1.0;
            for (unsigned int j = 0; j < dim; ++j) {
                if (i == j) {
                    A(i + dim, j) = -1.0;
                } else {
                    A(i + dim, j) = 0.0;
                }
            }
        }
    } else {

        m = 2 << (dim - 1);
        A.resize(m, dim);
        b.resize(m);
        for(unsigned int i=0; i<m; ++i){
            b(i) = 1;
            unsigned int k=i, j=0;
            while(k!=0){
                int bit = k % 2;
                if(bit == 0) {
                    A(i,j) = -1.0;
                }else {
                    A(i,j) = 1.0;
                }
                k = k >> 1;
                ++j;
            }
            for(; j<dim; ++j) {
                A(i,j) = -1.0;
            }
        }
    }
    Polytope P;
    P.init(dim, A, b);

    return P;
}


template <typename Polytope>
Polytope gen_cross(const unsigned int &dim, const bool &Vpoly) {

    unsigned int m;
    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;

    Polytope P;
    MT A;
    VT b;
    if (!Vpoly) {

        m = 2 << (dim - 1);
        A.resize(m, dim);
        b.resize(m);
        for(unsigned int i=0; i<m; ++i){
            b(i) = 1;
            unsigned int k=i, j=0;
            while(k!=0){
                unsigned int bit = k % 2;
                if(bit == 0) {
                    A(i,j) = -1.0;
                }else {
                    A(i,j) = 1.0;
                }
                k = k >> 1;
                ++j;
            }
            for(; j<dim; ++j) {
                A(i,j) = -1.0;
            }
        }
    } else {
        A.resize(2 * dim, dim);
        b.resize(2 * dim);

        for(unsigned int i=0; i<dim; ++i){
            b(i) = 1.0;
            for(unsigned int j=0; j<dim; ++j){
                if(i==j) {
                    A(i,j) = 1.0;
                } else {
                    A(i,j) = 0.0;
                }
            }
        }
        for(unsigned int i=0; i<dim; ++i){
            b(i + dim) = 1.0;
            for(unsigned int j=0; j<dim; ++j){
                if(i==j) {
                    A(i + dim, j) = -1.0;
                } else {
                    A(i + dim, j) = 0.0;
                }
            }
        }
    }
    P.init(dim, A, b);
    return P;
}


template <typename Polytope>
Polytope gen_simplex(const unsigned int &dim, const bool &Vpoly){
    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;

    MT A;
    A.resize(dim+1, dim);
    VT b;
    b.resize(dim+1);

    for(unsigned int i=0; i<dim; ++i){
        if (!Vpoly) {
            b(i) = 0;
        } else {
            b(i) = 1;
        }
        for(unsigned int j=0; j<dim; ++j){
            if(i==j) {
                A(i,j) = 1.0;
            } else {
                A(i,j) = 0.0;
            }
        }
    }
    b(dim) = 1;
    for(unsigned int j=0; j<dim; ++j){
        if (!Vpoly) {
            A(dim, j) = -1.0;
        } else {
            A(dim, j) = 0.0;
        }
    }
    Polytope P;
    P.init(dim, A, b);

    return P;
}


template <typename Polytope>
Polytope gen_prod_simplex(const unsigned int &dim, bool Vpoly = false){

    Polytope Perr;
    try
    {
        if(Vpoly) throw false;
    }
    catch (bool e) {
        #ifdef VOLESTI_DEBUG
        std::cout<<"Only prod simplices in H-representation can be generated.."<<std::endl;
        #endif
        return Perr;
    }

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;

    MT A;
    VT b;
    A.resize(2 * dim + 2, 2 * dim);
    b.resize(2 * dim + 2);
    Polytope P;

    //first simplex
    for(unsigned int i=0; i<dim; ++i){
        b(i) = 0.0;
        for(unsigned int j=0; j<dim; ++j) {
            A(i, j) = 0.0;
        }
        for(unsigned int j=0; j<dim; ++j){
            if(i==j) {
                A(i, j + dim) = 1.0;
            }else{
                A(i, j + dim) = 0.0;
            }
        }
    }

    b(dim) = 1.0;
    for(unsigned int j=0; j<dim; ++j) {
        A(dim, j) = 0.0;
    }
    for(unsigned int j=0; j<dim; ++j){
        A(dim, j + dim) = -1.0;
    }

    //second simplex
    for(unsigned int i=0; i<dim; ++i){
        b(dim + 1 + i) = 0.0;
        for(unsigned int j=0; j<dim; ++j){
            if(i==j) {
                A(dim + 1 + i, j) = 1.0;
            }else {
                A(dim + 1 + i, j) = 0.0;
            }
        }
        for(unsigned int j=0; j<dim; ++j) {
            A(dim + 1 + i, j + dim) = 0.0;
        }
    }
    b(2 * dim +1) = 1.0;
    for(unsigned int j=0; j<dim; ++j) {
        A(2 * dim +1, j) = -1.0;
    }
    for(unsigned int j=0; j<dim; ++j) {
        A(2 * dim +1, j + dim) = 0.0;
    }

    P.init(2 * dim, A, b);
    return P;
}


template <typename Polytope>
Polytope gen_skinny_cube(const unsigned int &dim, bool Vpoly = false) {

    Polytope Perr;
    try
    {
        if(Vpoly) throw false;
    }
    catch (bool e) {
        #ifdef VOLESTI_DEBUG
        std::cout<<"Only prod simplices in H-representation can be generated.."<<std::endl;
        #endif
        return Perr;
    }

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;

    MT A;
    A.resize(2 * dim, dim);
    VT b;
    b.resize(2 * dim);

    for(unsigned int i=0; i<dim; ++i){
        if (i==0) {
            b(i) = 100.0;
        } else {
            b(i) = 1.0;
        }
        for(unsigned int j=0; j<dim; ++j){
            if(i==j) {
                A(i,j) = 1.0;
            } else {
                A(i,j) = 0.0;
            }
        }
    }
    for(unsigned int i=0; i<dim; ++i){
        if (i==0) {
            b(i + dim) = 100.0;
        } else {
            b(i + dim) = 1.0;
        }
        for(unsigned int j=0; j<dim; ++j){
            if(i==j) {
                A(i + dim, j) = -1.0;
            } else {
                A(i + dim, j) = 0.0;
            }
        }
    }
    Polytope P;
    P.init(dim, A, b);

    return P;
}


/*
 * ToDo: brkhoff polytope generator
template <class Polytope>
Polytope gen_birk(int n, bool Vpoly = false){
 int m = pow(n,2);
  int d = pow(n-1,2)+1;

  std::cout << "birk_"<<n<<".ine\n";
  std::cout << "H-representation\n";
  std::cout << "begin\n";
  std::cout << " " << m << " " << d << " integer\n";

  std::cout << -1*(n-2) << " ";
  for(int j=1; j< d; ++j)
    std::cout << "1 ";
  std::cout << "\n";

  for(int i=0; i<n-1; ++i){
		std::cout << "1 ";
		for(int j=1; j< d; ++j){
			if(j%(n-1) == i){
		    std::cout << "-1 ";
		  }
		  else std::cout << " 0 ";
		}
		std::cout << "\n";
	}

	for(int i=0; i<n-1; ++i){
		std::cout << "1 ";
		for(int j=0; j< d-1; ++j){
			if(j/(n-1) == i){
		    std::cout << "-1 ";
		  }
		  else std::cout << " 0 ";
		}
		std::cout << "\n";
	}

	for(int i=0; i<d-1; ++i){
		std::cout << "0 ";
		for(int j=0; j< d-1; ++j){
			if(j == i){
		    std::cout << " 1 ";
		  }
		  else std::cout << " 0 ";
		}
		std::cout << "\n";
	}
	std::cout << "end\ninput_incidence" << std::endl;
 }


 */

#endif
