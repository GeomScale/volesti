// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef POLYTOPE_GENERATORS_H
#define POLYTOPE_GENERATORS_H

#include <exception>
#include "samplers.h"

template <class MT>
void removeRow(MT &matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.bottomRows(numRows-rowToRemove);

    matrix.conservativeResize(numRows,numCols);
}

template <class Polytope>
Polytope gen_cube(unsigned int dim, bool Vpoly) {

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


template <class Polytope>
Polytope gen_cross(unsigned int dim, bool Vpoly) {

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


template <class Polytope>
Polytope gen_simplex(unsigned int dim, bool Vpoly){
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


template <class Polytope>
Polytope gen_prod_simplex(unsigned int dim, bool Vpoly = false){

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


template <class Polytope>
Polytope gen_skinny_cube(unsigned int dim, bool Vpoly = false) {

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


template <class Polytope, class RNGType>
Polytope gen_zonotope(unsigned int dim, unsigned int m) {

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;
    typedef typename Polytope::NT    NT;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::normal_distribution<> rdist(0, 1);

    MT A;
    VT b;
    A.resize(m, dim);
    b.resize(m);
    Polytope P;

    for (unsigned int i = 0; i < m; ++i) {
        b(i) = 1.0;
        for (unsigned int j = 0; j < dim; ++j) {
            A(i,j) = rdist(rng);
        }
    }

    P.init(dim, A, b);
    return P;
}

template <class Polytope, class RNGType>
Polytope random_vpoly(unsigned int dim, unsigned int k) {

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;
    typedef typename Polytope::NT    NT;
    typedef typename Polytope::PolytopePoint Point;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);

    Point p;
    typename std::vector<NT>::iterator pit;
    MT V(k, dim);
    unsigned int j;

    for (unsigned int i = 0; i < k; ++i) {
        p = get_direction<RNGType, Point, NT>(dim);
        pit = p.iter_begin();
        j = 0;
        for ( ;  pit!=p.iter_end(); ++pit, ++j) {
            V(i,j) = *pit;
        }
    }

    Polytope VP;
    VT b = VT::Ones(k);
    VP.init(dim, V, b);

    return VP;

}


template <class Polytope, class RNGType>
Polytope random_vpoly_incube(unsigned int d, unsigned int k) {

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;
    typedef typename Polytope::NT    NT;
    typedef typename Polytope::PolytopePoint Point;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::random::uniform_real_distribution<> urdist1(-1, 1);

    Point p(d);
    typename std::vector<NT>::iterator pit;
    MT V(k, d);
    unsigned int j, count_row,it=0;
    std::vector<int> indices;
    Polytope VP;
    VT b = VT::Ones(k);

    for (unsigned int i = 0; i < k; ++i) {
        for (int j = 0; j < d; ++j) {
            V(i, j) = urdist1(rng);
        }
    }
    if(k==d+1){
        VP.init(d, V, b);
        return VP;
    }

    MT V2(k,d);
    V2 = V;
    indices.clear();
    while(it<20) {
        V.resize(V2.rows(), d);
        V = V2;
        for (int i = 0; i < indices.size(); ++i) {
            V.conservativeResize(V.rows()+1, d);
            for (int j = 0; j < d; ++j) {
                V(V.rows()-1, j) = urdist1(rng);
            }
        }
        indices.clear();
        V2.resize(k, d);
        V2 = V;

        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < d; ++j) {
                p.set_coord(j, V(i, j));
            }
            removeRow(V2, i);
            if (memLP_Vpoly(V2, p)){
                indices.push_back(i);
            }
            V2.resize(k, d);
            V2 = V;
        }
        if (indices.size()==0) {
            VP.init(d, V, b);
            return VP;
        }
        V2.resize(k - indices.size(), d);
        count_row =0;
        for (int i = 0; i < k; ++i) {
            if(std::find(indices.begin(), indices.end(), i) != indices.end()) {
                continue;
            } else {
                for (int j = 0; j < d; ++j) V2(count_row, j) = V(i,j);
                count_row++;
            }
        }
        it++;
    }

    VP.init(d, V2, VT::Ones(V2.rows()));
    return VP;

}


template <class Polytope, class RNGType>
Polytope random_hpoly(unsigned int dim, unsigned int m) {

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;
    typedef typename Polytope::NT    NT;
    typedef typename Polytope::PolytopePoint Point;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::random::uniform_real_distribution<> urdist1(-10, 10);
    Point p(dim);
    typename std::vector<NT>::iterator pit;
    MT A(m, dim);
    VT b(m);
    unsigned int j;

    for(unsigned int i=0; i<m; ++i){
        p = get_direction<RNGType, Point, NT>(dim);
        pit = p.iter_begin();
        j = 0;
        for ( ;  pit!=p.iter_end(); ++pit, ++j) {
            A(i,j) = *pit;
        }
        b(i) = 10.0;

    }
    Polytope HP;
    HP.init(dim, A, b);

    return HP;
}


template <class Polytope, class RNGType>
Polytope gen_zonotope_gaussian(unsigned int dim, unsigned int m) {

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;
    typedef typename Polytope::NT    NT;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::normal_distribution<> rdist(0, 1);
    boost::normal_distribution<> rdist2(50, 33.3);

    MT A;
    VT b;
    A.resize(m, dim);
    b.resize(m);
    Polytope P;
    NT rand_gaus;

    for (unsigned int i = 0; i < m; ++i) {
        b(i) = 1.0;
        for (unsigned int j = 0; j < dim; ++j) {
            A(i,j) = rdist(rng);
        }
        A.row(i)=A.row(i)/A.row(i).norm();
        while(true){
            rand_gaus = rdist2(rng);
            if (rand_gaus > 0.0 && rand_gaus<100.0){
                A.row(i) = A.row(i) * rand_gaus;
                break;
            }
        }
    }

    P.init(dim, A, b);
    return P;
}


template <class Polytope, class RNGType>
Polytope gen_zonotope_uniform(unsigned int dim, unsigned int m) {

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;
    typedef typename Polytope::NT    NT;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::normal_distribution<> rdist(0, 1);
    boost::random::uniform_real_distribution<> urdist1(0, 100);

    MT A;
    VT b;
    A.resize(m, dim);
    b.resize(m);
    Polytope P;

    for (unsigned int i = 0; i < m; ++i) {
        b(i) = 1.0;
        for (unsigned int j = 0; j < dim; ++j) {
            A(i,j) = rdist(rng);
        }
        A.row(i)=A.row(i)/A.row(i).norm();
        A.row(i) = A.row(i) * urdist1(rng);
    }

    P.init(dim, A, b);
    return P;
}


template <class Polytope, class RNGType>
Polytope gen_zonotope_exponential(unsigned int dim, unsigned int m) {

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;
    typedef typename Polytope::NT    NT;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::normal_distribution<> rdist(0, 1);
    boost::normal_distribution<> expdist(1.0/30.0);

    MT A;
    VT b;
    A.resize(m, dim);
    b.resize(m);
    Polytope P;
    NT rand_exp;

    for (unsigned int i = 0; i < m; ++i) {
        b(i) = 1.0;
        for (unsigned int j = 0; j < dim; ++j) {
            A(i,j) = rdist(rng);
        }
        A.row(i)=A.row(i)/A.row(i).norm();
        while(true){
            rand_exp = expdist(rng);
            if (rand_exp > 0.0 && rand_exp<100.0){
                A.row(i) = A.row(i) * rand_exp;
                break;
            }
        }
    }

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
