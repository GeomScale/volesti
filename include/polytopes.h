// RandGeom is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// RandGeom is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with HeaDDaCHe,
// see <http://www.gnu.org/licenses/>.
//
// Developer: Vissarion Fisikopoulos

#ifndef POLYTOPES_H
#define POLYTOPES_H

#define CGAL_QP_NO_ASSERTIONS

//this is for LP-solver
//#include "../external/kd_GeRaF/source/Auto_random_kd_forest.h"
#include "falconn/lsh_nn_table.h"
#include <iostream>
#include <chrono>
#include <Eigen/Eigen>
#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include <CGAL/intersections_d.h>
// choose exact integral type
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpzf.h>
#include <ANN/ANN.h>

typedef CGAL::Gmpzf ET;
#endif
//typedef double ET;
//#else
#include <CGAL/MP_Float.h>
//typedef CGAL::MP_Float ET;
//#endif
#include <CGAL/enum.h>
#include <boost/random/shuffle_order.hpp>
#include <rref.h>
//EXPERIMENTAL
//to implement boundary oracles using NN queries
//#include <flann/flann.hpp>

// my H-polytope class
template <typename K>
class stdHPolytope {
private:
    typedef std::vector<K>        stdCoeffs;
    typedef std::vector<stdCoeffs>  stdMatrix;
    //EXPERIMENTAL
    //typedef std::vector<flann::Index<flann::L2<double> > >  Flann_trees;

public:
	//stdHPolytope(const stdHPolytope&) e;
	stdHPolytope& operator=(const stdHPolytope&) = delete;
    stdHPolytope() {
		hptable = NULL;
		ANNkd_tree* kdTree = NULL;
		_maxDistToBoundary = 0;
	}

    // constructor: cube(d)
    stdHPolytope(int d): _d(d) {
        for(int i=0; i<d; ++i) {
            stdCoeffs coeffs;
            coeffs.push_back(K(1));
            for(int j=0; j<d; ++j) {
                if(i==j)
                    coeffs.push_back(K(1));
                else coeffs.push_back(K(0));
            }
            _A.push_back(coeffs);
        }
        for(int i=0; i<d; ++i) {
            stdCoeffs coeffs;
            coeffs.push_back(K(1));
            for(int j=0; j<d; ++j) {
                if(i==j)
                    coeffs.push_back(K(-1));
                else coeffs.push_back(K(0));
            }
            _A.push_back(coeffs);
        }
		hptable = NULL;
		ANNkd_tree* kdTree = NULL;
		_maxDistToBoundary = 0;
    }

	/**
	 * Polytope membership functions
	 */
	Point project(Point& point, int facet_idx) {
		/**
		 * proj=p - a * ((p.dot_product(a) + b)/ a.dot_product(a))
		 */
		typename stdCoeffs::iterator mit=_A[facet_idx].begin();
		++mit;
		Vector facet_normal = Vector(_d, mit, _A[facet_idx].end());
		Vector point_as_v = point - CGAL::ORIGIN;

		double dir_coeff = ((point_as_v * facet_normal) + _A[facet_idx][0]) / (facet_normal * facet_normal);
		Vector facet_normal_mul(facet_normal);
		facet_normal_mul *= dir_coeff;
		Vector projection = point_as_v - facet_normal_mul;

		return CGAL::ORIGIN + projection;
	}

	Point _get_reflexive_point(Point& internalPoint, int facet_idx) {
		Point projection = this->project(internalPoint, facet_idx);
		Vector projection_as_v = projection - CGAL::ORIGIN;
		Vector internalPoint_as_v = internalPoint - CGAL::ORIGIN;

		//(2 * projection - internalPoint); // <=> internalPoint + 2*(projection-internalPoint)
		projection_as_v *= 2;
		Vector reflexive_point = projection_as_v - internalPoint_as_v;

		return CGAL::ORIGIN + reflexive_point;
	}

	void center_sites(std::vector<Point>& sites, Point& internalPoint) {
		Vector internalPoint_as_v = internalPoint - CGAL::ORIGIN;
		for (int i=0; i<sites.size(); i++) {
			_sites.push_back(
				sites[i]
				//(CGAL::ORIGIN + ((sites[i] - CGAL::ORIGIN) - internalPoint_as_v))
			);
		}
		_sites.push_back(CGAL::ORIGIN + (internalPoint_as_v - internalPoint_as_v));
	}

	Point create_point_representation() {
		data.clear();
		_sites.clear();
		Point tmp_internalPoint;
		double radius;

		this->chebyshev_center(tmp_internalPoint, radius);
		auto it = tmp_internalPoint.cartesian_begin();
		++it;
		Point internalPoint(_d, it, tmp_internalPoint.cartesian_end());
		std::cout << "Chebyshev center: " << internalPoint << std::endl;
		std::vector<Point> sites;

		for (int i=0; i<this->_A.size(); i++) {
			Point reflexive_point = this->_get_reflexive_point(internalPoint, i);
			sites.push_back(reflexive_point);
			double tmpDist = ((reflexive_point-CGAL::ORIGIN)-(internalPoint-CGAL::ORIGIN)).squared_length(); 
			if (tmpDist > _maxDistToBoundary) {
				_maxDistToBoundary = tmpDist;
			}
		}
		center_sites(sites, internalPoint); // This is useful for the FALCON LSH structure

		sites.clear();
		return _sites[_sites.size()-1];
	}

	void create_ann_jl_ds() {
		int new_d = (int)std::ceil(std::log2(_sites.size())/(std::log2(std::log2(_sites.size()))));
		std::default_random_engine generator;
		std::normal_distribution<double> distribution(0.0,1.0);
		proj_matrix.resize(dimension(), new_d);
		for (int i=0; i<dimension(); i++) {
			for (int j=0; j<new_d; j++) {
				proj_matrix(i, j) = distribution(generator);
			}
		}
		Eigen::MatrixXd map(_sites.size(), dimension());
		for (int i = 0; i < _sites.size(); i++) {
			auto it = _sites[i].cartesian_begin();	
			for (int j=0; j<dimension(); j++, ++it) {
				map(i,j) = (*it);
			}
		}
		Eigen::MatrixXd projected_data = map*proj_matrix;

		ANNpointArray point_array = annAllocPts(_sites.size(), new_d);
		for (int i=0; i<_sites.size(); i++) {
			for (int j=0; j<new_d; j++) {
				point_array[i][j] = projected_data(i,j);
			}
		}
		this->bdTree = new ANNbd_tree(
						point_array,
						_sites.size(),
						new_d
		);
	}

	void create_lsh_ds(int k, int l) {

		falconn::LSHConstructionParameters params_hp;
		this->_k = k;
		this->_l = l;

		uint64_t seed = 119417657;
		params_hp.dimension = this->dimension();
		params_hp.lsh_family = falconn::LSHFamily::Hyperplane;
		params_hp.distance_function = falconn::DistanceFunction::EuclideanSquared;
		params_hp.storage_hash_table = falconn::StorageHashTable::FlatHashTable;
		params_hp.k = k;
		params_hp.l = l;
		params_hp.num_setup_threads = 0;
		params_hp.seed = seed ^ 833840234;

		for (int i=0; i<_sites.size(); i++) {
			std::vector<double> tmp_vec(_sites[i].cartesian_begin(), _sites[i].cartesian_end());
			Eigen::Map<Eigen::VectorXd> map(&tmp_vec[0], this->_d);
			this->data.push_back(map);
		}

		this->hptable = falconn::construct_table<falconn::DenseVector<double>>(this->data, params_hp);
	}

	void create_ann_ds() {
		ANNpointArray point_array = annAllocPts(_sites.size(), dimension());
		for (int i=0; i<_sites.size(); i++) {
			for (int j=0; j<dimension(); j++) {
				point_array[i][j] = _sites[i][j];
			}
		}
		this->kdTree = new ANNkd_tree(
						point_array,
						_sites.size(),
						this->dimension());
	}

	bool contains_point_ann_jl(ANNpoint& p, ANNidxArray& nnIdx, ANNdistArray& dists, double membership_epsilon, int* nnIndex_ptr) {
		int new_d = (int)std::ceil(std::log2(_sites.size())/(std::log2(std::log2(_sites.size()))));
		auto queryPt = annAllocPt(new_d);
		for (int i=0; i<new_d; i++) {
			for (int j=0; j<dimension(); j++) {
				queryPt[i] += p[j] * proj_matrix(j, i);
			}
		}

		int size = (int)std::ceil(std::pow(1+num_of_hyperplanes(), (double)1/4));
		//int size = (1+num_of_hyperplanes())/100;//i(int)std::ceil(std::pow(1+num_of_hyperplanes(), (double)1/3));

		
		double epsilon = (2*membership_epsilon)/(1-membership_epsilon);
		this->bdTree->annkSearch(
				queryPt,
				size,
				nnIdx,
				dists,
				epsilon
		);

		(*nnIndex_ptr) = -1;
		double minDist = 1000000;
		bool contains_center = false;
		for (int i=0; i<size; i++) {
			double sum = 0;
			auto it = _sites[nnIdx[i]].cartesian_begin();
			for (int j=0; j<dimension(); j++, ++it) {
				sum += std::pow(p[j] - (*it), 2);
			}
			if (sum<minDist) {
				minDist = sum;
				(*nnIndex_ptr) = nnIdx[i];
			}
			if (nnIdx[i]==_sites.size()-1) {
				contains_center = true;
			}
		}

		if ((*nnIndex_ptr)!=_sites.size()-1) {
			double sum = 0;
			auto it = _sites[_sites.size()-1].cartesian_begin();
			for (int j=0; j<dimension(); j++, ++it) {
				sum += std::pow(p[j] - (*it), 2);
			}
			//std::cout << "center dist " << sum << " -- ann dist " << minDist << std::endl;
			if (sum<minDist) {
				(*nnIndex_ptr) = _sites.size()-1;
			}
		}

		return (*nnIndex_ptr)==_sites.size()-1;
	}

	bool contains_point_ann(ANNpoint p, ANNidxArray nnIdx, ANNdistArray dists, double membership_epsilon, int* nnIndex_ptr) {
		//auto it = _sites[_sites.size()-1].cartesian_begin();
		//for (int i=0; i<dimension(); i++, ++it) {
		//	p[i] -= (*it);
		//}	
		double epsilon = (2*membership_epsilon)/(1-membership_epsilon);
		this->kdTree->annkSearch(
				p,
				1,
				nnIdx,
				dists,
				epsilon
		);

		bool is_in = nnIdx[0]==_sites.size()-1;
		(*nnIndex_ptr) = nnIdx[0];
		return is_in;
	}

	bool contains_point_lsh(Point& p, int num_probes, int* nnIndex_ptr) {
		Point newP = (CGAL::ORIGIN + (p-CGAL::ORIGIN)-(_sites[_sites.size()-1]-CGAL::ORIGIN));
		this->hptable->set_num_probes(num_probes);
		std::vector<double> tmp_vec(newP.cartesian_begin(), newP.cartesian_end());
		Eigen::Map<Eigen::VectorXd> map(&tmp_vec[0], this->_d);
		(*nnIndex_ptr) = this->hptable->find_nearest_neighbor(map);
		//std::cout << "nn index: " << (*nnIndex_ptr) << std::endl;
		return  (*nnIndex_ptr)== _sites.size()-1;
	}

    bool contains_point_naive(Point p, double epsilon) {
        //std::cout << "Running is in" << std::endl;
        //exit(1);
		int idx = 0;
        for(typename stdMatrix::iterator mit=_A.begin(); mit<_A.end(); ++mit) {
            typename stdCoeffs::iterator lit;
            Point::Cartesian_const_iterator pit;
            pit=p.cartesian_begin();
            lit=mit->begin();
			K coeff = (*lit);
            K sum=0;//(*lit);
            ++lit;
            for( ; lit<mit->end() ; ++lit, ++pit) {
                //std::cout << *lit << " " << *pit <<std::endl;
                sum += *lit * (*pit);
            }

            //std::cout<<sum<<std::endl;
            if (sum>coeff-epsilon) {
				//auto lit2 = mit->begin();
				//coeff = (*lit2);
				//lit2++;
				//std::cout << "SUM: " << sum << " (";
				//for (; lit2!=mit->end(); lit2++) {
				//	std::cout << (*lit2) << ", ";
				//}
				//std::cout  << ") <= " << coeff << "  is violated" << std::endl;
                return false;
			}
        }
        return true;
    }
	/**
	 * End of polytope membership functions
	 */

	/**
	 * Polytope boundary functions
	 */
	Point compute_boundary_intersection(Point& point, Vector& vector, double epsilon, bool use_lsh) {
		return compute_boundary_intersection(Ray(point, vector), epsilon, use_lsh);
	}

	Point compute_boundary_intersection(Ray& ray, double epsilon, bool use_jl) {
		Line ray_line(ray);
		Vector ray_direction = ray.direction().vector();		
		ray_direction *= dimension() * _maxDistToBoundary;
		Point x0 = CGAL::ORIGIN + ((ray.source()-CGAL::ORIGIN) + ray_direction);
		int nnIndex = -1;
		//std::cout << "x0 is in ? " << contains_point_naive(x0, 0) << std::endl;

		Vector ray_source_v = (ray.source()-CGAL::ORIGIN);//- (_sites[_sites.size()-1]-CGAL::ORIGIN);
		ANNidxArray annIdx;
		ANNdistArray dists;

		int size = (int)std::ceil(std::pow(1+num_of_hyperplanes(), (double)1/4));
		annIdx = new ANNidx[size];
		dists = new ANNdist[size];
		ANNpoint queryPt = annAllocPt(dimension());

		do {
			auto start_time = std::chrono::high_resolution_clock::now(); 
			bool contains;
			auto it = x0.cartesian_begin();
			for (int i=0; it<x0.cartesian_end(); i++, it++) {
				queryPt[i] = (*it);
			}
			if (!use_jl) {
				contains = contains_point_ann(queryPt, annIdx, dists, epsilon, &nnIndex);
			} else {
				contains = contains_point_ann_jl(queryPt, annIdx, dists, epsilon, &nnIndex);
			}
			auto end_time = std::chrono::high_resolution_clock::now();
			auto elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
			double elapsed_total = elapsed.count();
			//std::cout << "ANN took: " << elapsed.count() << "s" << std::endl;
			//bool naive_contains = contains_point_naive(x0, 0);
			//int i=0;
			//Point projection = project(x0, i++);
			//double minDist = std::sqrt(((x0-CGAL::ORIGIN)-(projection-CGAL::ORIGIN)).squared_length());
			//for (; i<num_of_hyperplanes(); i++) {
			//	projection = project(x0, i);
			//	double tmpDist = std::sqrt(((x0-CGAL::ORIGIN)-(projection-CGAL::ORIGIN)).squared_length());
			//	if (tmpDist<minDist) {
			//		minDist = tmpDist;
			//	}
			//}
			//std::cout << "contained: " << contains<< " -- Contained naive: " << naive_contains << " minDist: " << minDist <<" -- LSH NN: " << nnIndex << " -- Point: " << x0 << std::endl;
		
			if (!contains) {
				start_time = std::chrono::high_resolution_clock::now(); 

				auto it = _A[nnIndex].begin();
				double coeff = (*it);
				++it;
				Hyperplane nn_facet(dimension(), it, _A[nnIndex].end(), coeff);
				CGAL::cpp11::result_of<Kernel::Intersect_d(Line, Hyperplane)>::type x1_tmp = CGAL::intersection(ray_line, nn_facet);
				Point* x1 = boost::get<Point>(&*x1_tmp);
				double x1_ray_norm = (((*x1)-CGAL::ORIGIN) - (ray_source_v)).squared_length();
				double x0_ray_norm = ((x0-CGAL::ORIGIN) - (ray_source_v)).squared_length();
				end_time = std::chrono::high_resolution_clock::now();
				elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
				elapsed_total += elapsed.count();
				//std::cout << "Intersection took: " << elapsed.count() << "s" << std::endl;
				//std::cout << "On ray: " << ray.has_on((*x1)) << " -- " << x0_ray_norm << "\t" << x1_ray_norm << std::endl;
				if (!ray.has_on((*x1)) || x1_ray_norm>=x0_ray_norm) {
					start_time = std::chrono::high_resolution_clock::now(); 
					Vector newPoint_v = ((x0-CGAL::ORIGIN) - (ray.source()-CGAL::ORIGIN));
					newPoint_v *= (1-0.5);//*_maxDistToBoundary);
					(*x1) = CGAL::ORIGIN + (newPoint_v + ray_source_v);
					//std::cout << "\tNew norm: " << (((*x1)-CGAL::ORIGIN)-ray_source_v).squared_length() << std::endl;
					end_time = std::chrono::high_resolution_clock::now();
					elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
					elapsed_total += elapsed.count();
					//std::cout << "New x1 took: " << elapsed.count() << "s" << std::endl;
				}
				x0 = Point(dimension(), (*x1).cartesian_begin(), (*x1).cartesian_end());
			}
			//std::cout << "1 iteration took (approximately) " << elapsed_total << "s" << std::endl;
			
		} while (nnIndex!=_sites.size()-1);
		delete []queryPt;
		delete []annIdx;
		delete []dists;

		return x0;
	}
	/**
	 * End of polytope boundary functions
	 */


    int dimension() {
        return _d;
    }

    int num_of_hyperplanes() {
        return _A.size();
    }

    K get_coeff(int i, int j) {
        return _A[i][j];
    }

    void put_coeff(int i, int j, K value) {
        _A[i][j] = value;
    }

    // default initialize: cube(d)
    int init(int d) {
        _d=d;
        for(int i=0; i<d; ++i) {
            stdCoeffs coeffs;
            coeffs.push_back(K(1));
            for(int j=0; j<d; ++j) {
                if(i==j)
                    coeffs.push_back(K(1));
                else coeffs.push_back(K(0));
            }
            _A.push_back(coeffs);
        }
        for(int i=0; i<d; ++i) {
            stdCoeffs coeffs;
            coeffs.push_back(K(1));
            for(int j=0; j<d; ++j) {
                if(i==j)
                    coeffs.push_back(K(-1));
                else coeffs.push_back(K(0));
            }
            _A.push_back(coeffs);
        }
        return 0;
    }

    int init(stdMatrix Pin) {
        _d = Pin[0][1]-1;
        typename stdMatrix::iterator pit=Pin.begin();
        ++pit;
        for( ; pit<Pin.end(); ++pit) {
            _A.push_back(*pit);
        }
        //double tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        //std::random_shuffle (_A.begin(), _A.end());
        //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        //std::shuffle (_A.begin(), _A.end(), std::default_random_engine(seed));
        //std::shuffle (_A.begin(), _A.end(), std::default_random_engine(seed));
        //boost::random::shuffle_order_engine<stdMatrix>(_A.begin(), _A.end());
        //double tstop = (double)clock()/(double)CLOCKS_PER_SEC;
        //std::cout << "Shuffle time = " << tstop - tstart << std::endl;
        return 0;
    }

    // print polytope in input format
    int print() {
        std::cout<<" "<<_A.size()<<" "<<_d+1<<" float"<<std::endl;
        for(typename stdMatrix::iterator mit=_A.begin(); mit<_A.end(); ++mit) {
            for(typename stdCoeffs::iterator lit=mit->begin(); lit<mit->end() ; ++lit)
                std::cout<<*lit<<" ";
            std::cout<<std::endl;
        }
        return 0;
    }

    /*
    int is_in(stdCoeffs p) {
    	for(typename stdMatrix::iterator mit=_A.begin(); mit<_A.end(); ++mit){
    		typename stdCoeffs::iterator lit,pit;
    		pit=p.begin();
    		lit=mit->coefficients_begin();
    		K sum=(*lit);
    		++lit;
    		for( ; lit<mit->coefficients_end() ; ++lit){
    			//std::cout << *lit << " " << *pit <<std::endl;
    			sum += *lit * (*pit);
    		}
    		//std::cout<<sum<<std::endl;
    		if(sum<K(0))
    			return mit-_A.begin();
    	}
    	return -1;
    }
    */

    // Compute the reduced row echelon form
    // used to transofm {Ax=b,x>=0} to {A'x'<=b'}
    // e.g. Birkhoff polytopes
    int rref() {
        to_reduced_row_echelon_form(_A);
        std::vector<int> zeros(_d+1,0);
        std::vector<int> ones(_d+1,0);
        std::vector<int> zerorow(_A.size(),0);
        for (int i = 0; i < _A.size(); ++i)
        {
            for (int j = 0; j < _d+1; ++j) {
                if ( _A[i][j] == double(0)) {
                    ++zeros[j];
                    ++zerorow[i];
                }
                if ( _A[i][j] == double(1)) {
                    ++ones[j];
                }
            }
        }
        for(typename stdMatrix::iterator mit=_A.begin(); mit<_A.end(); ++mit) {
            int j =0;
            for(typename stdCoeffs::iterator lit=mit->begin(); lit<mit->end() ; ) {
                if(zeros[j]==_A.size()-1 && ones[j]==1)
                    (*mit).erase(lit);
                else { //reverse sign in all but the first column
                    if(lit!=mit->end()-1) *lit = (-1)*(*lit);
                    ++lit;
                }
                ++j;
            }
        }
        //swap last and first columns
        for(typename stdMatrix::iterator mit=_A.begin(); mit<_A.end(); ++mit) {
            double temp=*(mit->begin());
            *(mit->begin())=*(mit->end()-1);
            *(mit->end()-1)=temp;
        }
        //delete zero rows
        for (typename stdMatrix::iterator mit=_A.begin(); mit<_A.end(); ) {
            int zero=0;
            for(typename stdCoeffs::iterator lit=mit->begin(); lit<mit->end() ; ++lit) {
                if(*lit==double(0)) ++zero;
            }
            if(zero==(*mit).size())
                _A.erase(mit);
            else
                ++mit;
        }
        //update _d
        _d=(_A[0]).size();
        // add unit vectors
        for(int i=1; i<_d; ++i) {
            std::vector<double> e(_d,0);
            e[i]=1;
            _A.push_back(e);
        }
        // _d should equals the dimension
        _d=_d-1;
        return 1;
    }

    /* EXPERIMENTAL

    // compute the dual representation of P
    // and construct the flann trees used for NN queries
    int dual(int dir){
    	int d=_d;
    	//std::cout<<"\nDual\n";
    	for (int k=1; k<=d; ++k){
    	//for (int k=d; k>=1; --k){
    		flann::Matrix<double> dataset(new double[_A.size()*(d+1)], _A.size(), (d+1));
    		stdMatrix _Adual;
    		double M=0;
    		for (int i = 0; i < _A.size(); ++i)
    		{
    			double t_norm_squared=0;
    			stdCoeffs coeffs;
    			for (int j = 0; j < d+1; ++j){
    	      if (j != k){
    					double value = dir*_A[i][j]/(_A[i][k]);
    					coeffs.push_back(value);
    					//datasets[k-1][i][j] = value;
    					t_norm_squared += std::pow(value,2);
    					//std::cout<<value<<" ";
    				}
    	    }
    	    for (int j = 0; j < d-1; ++j){
    				dataset[i][j] = coeffs[j+1];
    			}
    			dataset[i][d-1]=-1*coeffs[0];
    	    //coeffs.push_back(t_norm_squared);
    	    dataset[i][d]=t_norm_squared;
    	    if(1+t_norm_squared > M) M=1+t_norm_squared;
    	    //_Adual.push_back(coeffs);
    			//std::cout<<"\n";
    	  }
    	  for (int i = 0; i < _A.size(); ++i)
    		{
    			//_Adual[i].push_back(std::sqrt(M-_Adual[i][d+1]));
    			dataset[i][d] = std::sqrt(M-1-dataset[i][d]);
    			//std::cout<<M-1-_Adual[i][d+1]<<" ";
    		}
    	  //std::cout<<"\n-----\n";
    		//_Aduals.push_back(_Adual);

    		//print dataset k

    		//for (int i = 0; i < _A.size(); ++i){
    		//  for (int j = 0; j < d+1; ++j)
    		//		std::cout<<dataset[i][j]<<" ";
    		//  std::cout<<"\n";
    		//}
    		// construct an randomized kd-tree index using  kd-trees
    		//flann::Index<flann::L2<double> > index(dataset, flann::KDTreeSingleIndexParams(4));
    		flann::Index<flann::L2<double> > index(dataset, flann::LinearIndexParams());
    		//Inexact results & Seg.faults
    		//flann::Index<flann::L2<double> > index(dataset, flann::KDTreeIndexParams(4));
    		index.buildIndex();
    		flann_trees.push_back(index);
    		//std::cout<<"\n-----\n";
    	}
    	return 1;
    }



    // dual must be set before calling this
    std::pair<NT,NT>
    query_dual(Point p, int rand_coord){
    	std::pair<NT,NT> NNpair;
    	int checks = 16;
    	int eps = 0;
    	int d = _d;
    	int Q=1;
    	int k=1;
    	flann::Matrix<int> indices(new int[Q*k], Q, k);
    	flann::Matrix<double> query(new double[Q*d], Q, d);
    	flann::Matrix<double> dists(new double[Q*k], Q, k);

    	// transform point for query
    	int j=0;
    	for(int i=0;i<d;i++){
    		if(i!=rand_coord)
    		  query[0][j++]=p.cartesian(i);
    	}
    	query[0][d-1]=-1;
    	query[0][d]=0;

    	//farthest
    	{
    	flann_trees[rand_coord].knnSearch(query, indices, dists, k,
    	                                  flann::SearchParams(checks, eps));

    	//std::cout<<"farthest rand_coord: "<<rand_coord<<std::endl;
    	//std::cout<<"Qresult: "<<indices[0][0]<<std::endl;
    	//std::cout<<"Qdist: "<<dists[0][0]<<std::endl;

    	int i = indices[0][0];
    	//
    	K lamda;
    	typename stdCoeffs::iterator cit;
    	Point::Cartesian_const_iterator rit;
    	rit=p.cartesian_begin();
    	cit=_A[i].begin();
    	K sum_nom=(*cit);
    	++cit;
    	K sum_denom= *(cit+rand_coord);
    	//std::cout<<ait->begin()-ait->end()<<" "<<r.cartesian_begin()-r.cartesian_end()<<" "<<
    	//         std::endl;
    	for( ; cit < _A[i].end() ; ++cit, ++rit){
    		sum_nom -= *cit * (*rit);
    	}
    	//lamdas[ait-_A.begin()] = sum_nom;
    	if(sum_denom==K(0)){
        //std::cout<<"div0"<<std::endl;
        ;
      }
      else{
        lamda = sum_nom*(1/sum_denom);
    	}
    	//std::cout<<"Qlambda:"<<lamda<<std::endl;
    	NNpair.first = lamda;
      }

    	//nearest
    	{
    	flann_trees[rand_coord+d].knnSearch(query, indices, dists, k,
    	                                  flann::SearchParams(checks, eps));

    	//std::cout<<"nearest rand_coord: "<<rand_coord<<std::endl;
    	//std::cout<<"Qresult: "<<indices[0][0]<<std::endl;
    	//std::cout<<"Qdist: "<<dists[0][0]<<std::endl;

    	int i = indices[0][0];
    	//
    	K lamda;
    	typename stdCoeffs::iterator cit;
    	Point::Cartesian_const_iterator rit;
    	rit=p.cartesian_begin();
    	cit=_A[i].begin();
    	K sum_nom=(*cit);
    	++cit;
    	K sum_denom= *(cit+rand_coord);
    	//std::cout<<ait->begin()-ait->end()<<" "<<r.cartesian_begin()-r.cartesian_end()<<" "<<
    	//         std::endl;
    	for( ; cit < _A[i].end() ; ++cit, ++rit){
    		sum_nom -= *cit * (*rit);
    	}
    	//lamdas[ait-_A.begin()] = sum_nom;
    	if(sum_denom==K(0)){
        //std::cout<<"div0"<<std::endl;
        ;
      }
      else{
        lamda = sum_nom*(1/sum_denom);
    	}
    	//std::cout<<"Qlambda:"<<lamda<<std::endl;
      NNpair.second = lamda;
      }

      // deallocate memory
      //delete[] query.ptr();
    //delete[] indices.ptr();
    //delete[] dists.ptr();
      //exit(1);
    	return NNpair;
    }
    */

    int is_in(Point p) {
        //std::cout << "Running is in" << std::endl;
        //exit(1);
        for(typename stdMatrix::iterator mit=_A.begin(); mit<_A.end(); ++mit) {
            typename stdCoeffs::iterator lit;
            Point::Cartesian_const_iterator pit;
            pit=p.cartesian_begin();
            lit=mit->begin();
            K sum=(*lit);
            ++lit;
            for( ; lit<mit->end() ; ++lit, ++pit) {
                //std::cout << *lit << " " << *pit <<std::endl;
                sum += *lit * (*pit);
            }

            //std::cout<<sum<<std::endl;
            if(sum<K(0))
                return mit-_A.begin();
        }
        return -1;
    }

    int chebyshev_center(Point& center, double& radius) {
        typedef CGAL::Linear_program_from_iterators
        <K**,                             // for A
        K*,                              // for b
        CGAL::Const_oneset_iterator<CGAL::Comparison_result>,  // for r
        bool*,                           // for fl
        K*,                              // for l
        bool*,                           // for fu
        K*,                              // for u
        K*>                              // for c
        Program;
        typedef CGAL::Quadratic_program_solution<ET> Solution;

        //std::cout<<"Cheb"<<std::endl;
        K* b = new K[_A.size()];
        K** A_col = new K*[_d+1];
        for(size_t i = 0; i < _d+1; ++i)
            A_col[i] = new K[_A.size()];

        stdMatrix B(_A);
        std::random_shuffle (B.begin(), B.end());
        for(size_t i=0; i<B.size(); ++i) {
            K sum_a2 = 0;
            b[i] = B[i][0];
            for(size_t j=0; j<_d; ++j) {
                A_col[j][i] = B[i][j+1];
                sum_a2 += std::pow(B[i][j+1],2);
            }
            A_col[_d][i] = std::sqrt(sum_a2);
        }

        CGAL::Const_oneset_iterator<CGAL::Comparison_result>
        r(CGAL::SMALLER);

        bool* fl = new bool[_d+1]();
        bool* fu = new bool[_d+1]();

        for(size_t i=0; i<_d+1; ++i)
            fl[i]=false;
        for(size_t i=0; i<_d+1; ++i)
            fu[i]=false;
        fl[_d]=true;

        K* l = new K[_d+1]();
        K* u = new K[_d+1]();
        K* c = new K[_d+1]();

        for(size_t i=0; i<_d+1; ++i)
            l[i]=K(0);
        for(size_t i=0; i<_d+1; ++i)
            u[i]=K(0);
        for(size_t i=0; i<_d+1; ++i)
            c[i]=K(0);
        c[_d]=K(-1);

        Program lp (_d+1, int(_A.size()), A_col, b, r,
                    fl, l, fu, u, c, 0);

        //CGAL::Quadratic_program_options options;
        //options.set_verbosity(1);                         // verbose mode
        //options.set_pricing_strategy(CGAL::QP_BLAND);     // Bland's rule
        //options.set_auto_validation(true);                // automatic self-check
        //Solution s = CGAL::solve_linear_program(lp, ET(), options);

        Solution s = CGAL::solve_linear_program(lp, ET());
        // output solution
        //std::cout << s;
        if (s.is_infeasible()) {
            std::cout << "The polytope P is unbounded and Vol(P)=0\n";
            exit(-1);
        }
        else {
            assert (s.is_optimal());
            Solution::Variable_value_iterator it = s.variable_values_begin();
            std::vector<double> vecp;
            for(; it!=s.variable_values_end()-1; ++it) {
                //std::cout<<CGAL::to_double(*it)<<" ";
                vecp.push_back(CGAL::to_double(*it));
            }
            center = Point(_d,vecp.begin(),vecp.end());
            //std::cout << center;
            radius = CGAL::to_double(*it);
            //std::cout << radius << std::endl;

        }
        // deallocate memory
        delete [] l;
        delete [] u;
        delete [] c;
        delete [] fl;
        delete [] fu;
        for(size_t i = 0; i < _d+1; ++i)
            delete [] A_col[i];
        delete [] b;
        return 0;
    }

    // compute intersection point of ray starting from r and pointing to v
    // with polytope discribed by _A
    std::pair<Point,Point> line_intersect(Point r,
                                          Vector v) {
        //std::cout<<"line-polytope-intersection"<<std::endl;
        K lamda=0;
        K min_plus=0, max_minus=0;
        bool min_plus_not_set=true;
        bool max_minus_not_set=true;
        for(typename stdMatrix::iterator ait=_A.begin(); ait<_A.end(); ++ait) {
            typename stdCoeffs::iterator cit;
            Point::Cartesian_const_iterator rit;
            rit=r.cartesian_begin();
            Point::Cartesian_const_iterator vit;
            vit=v.cartesian_begin();
            cit=ait->begin();
            K sum_nom=(*cit);
            ++cit;
            K sum_denom=K(0);
            //std::cout<<ait->begin()-ait->end()<<" "<<r.cartesian_begin()-r.cartesian_end()<<" "<<
            //         v.cartesian_begin()-v.cartesian_end()<<std::endl;
            for( ; cit < ait->end() ; ++cit, ++rit, ++vit) {
                //std::cout << sum_nom << " " << sum_denom <<std::endl;
                //std::cout << int(rit-r.cartesian_begin()) << " " << int(vit-v.cartesian_begin()) <<std::endl;
                sum_nom -= *cit * (*rit);
                sum_denom += *cit * (*vit);
            }
            if(sum_denom==K(0)) {
                std::cout<<"div0"<<std::endl;
            }
            else {
                lamda = sum_nom/sum_denom;
                if(min_plus_not_set && lamda>0) {
                    min_plus=lamda;
                    min_plus_not_set=false;
                }
                if(max_minus_not_set && lamda<0) {
                    max_minus=lamda;
                    max_minus_not_set=false;
                }
                if(lamda<min_plus && lamda>0) min_plus=lamda;
                if(lamda>max_minus && lamda<0) max_minus=lamda;
            }
            //std::cout<<r+(lamda*v)<<"\n"<<lamda<<std::endl;
        }
        /*
        std::cout<<"lmin,lmax= "<<max_minus<<" "<<min_plus<<std::endl;
        std::cout<<"r= "<<r<<std::endl;
        std::cout<<"v= "<<v<<std::endl;
        std::cout<<"p1= "<<r+(min_plus*v)<<std::endl;
        std::cout<<"p2= "<<r+(max_minus*v)<<std::endl;
        */
        return std::pair<Point,Point> (r+(min_plus*v),r+(max_minus*v));
    }

    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          int rand_coord) {
        //std::cout<<"line-polytope-intersection"<<std::endl;
        K lamda=0;
        //std::vector<NT> new_lamdas(_A.size());
        //std::vector<NT> new_lamdas;
        K min_plus=0, max_minus=0;
        bool min_plus_not_set=true;
        bool max_minus_not_set=true;
        //std::vector<NT>::iterator lamdait = lamdas.begin();
        for(typename stdMatrix::iterator ait=_A.begin(); ait<_A.end(); ++ait) {
            typename stdCoeffs::iterator cit;
            Point::Cartesian_const_iterator rit;
            rit=r.cartesian_begin();
            //Point::Cartesian_const_iterator vit;
            //vit=v.cartesian_begin();
            cit=ait->begin();
            K sum_nom=(*cit);
            ++cit;
            //here we just need the "rand_coord" coordinate of c
            //std::cout<<*(cit+rand_coord)<<"c= "<<std::endl;
            //for(typename stdCoeffs::iterator cit2=ait->begin() ; cit2 < ait->end() ; ++cit2){
            //  std::cout<<*cit2<<" ";
            //}
            //std::cout<<std::endl;
            K sum_denom= *(cit+rand_coord);
            //std::cout<<ait->begin()-ait->end()<<" "<<r.cartesian_begin()-r.cartesian_end()<<" "<<
            //         v.cartesian_begin()-v.cartesian_end()<<std::endl;
            for( ; cit < ait->end() ; ++cit, ++rit) {
                //std::cout << sum_nom << " " << sum_denom <<std::endl;
                //std::cout << int(rit-r.cartesian_begin()) << " " << int(vit-v.cartesian_begin()) <<std::endl;
                sum_nom -= *cit * (*rit);
                //sum_denom += *cit * (*vit);
            }
            //std::cout << sum_nom << " / "<< sum_denom<<std::endl;
            if(sum_denom==K(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            }
            else {
                lamda = sum_nom*(1/sum_denom);
                //lamdas[ait-_A.begin()] = lamda;

                if(min_plus_not_set && lamda>0) {
                    min_plus=lamda;
                    min_plus_not_set=false;
                }
                if(max_minus_not_set && lamda<0) {
                    max_minus=lamda;
                    max_minus_not_set=false;
                }
                if(lamda<min_plus && lamda>0) min_plus=lamda;
                if(lamda>max_minus && lamda<0) max_minus=lamda;
            }
            //std::cout<<r+(lamda*v)<<"\n"<<lamda<<std::endl;
        }
        return std::pair<NT,NT> (min_plus,max_minus);
    }

    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          Point &r_prev,
                                          int rand_coord,
                                          int rand_coord_prev,
                                          std::vector<NT> &lamdas,
                                          bool init) {
        //std::cout<<"line-polytope-intersection"<<std::endl;
        K lamda=0;
        std::vector<NT>::iterator lamdait = lamdas.begin();

        K min_plus=0, max_minus=0;
        bool min_plus_not_set=true;
        bool max_minus_not_set=true;
        int mini, maxi;

        if(init) { //first time compute the innerprod cit*rit
            for(typename stdMatrix::iterator ait=_A.begin(); ait<_A.end(); ++ait) {
                typename stdCoeffs::iterator cit;
                Point::Cartesian_const_iterator rit;
                rit=r.cartesian_begin();
                cit=ait->begin();
                K sum_nom=(*cit);
                ++cit;
                K sum_denom= *(cit+rand_coord);
                //std::cout<<ait->begin()-ait->end()<<" "<<r.cartesian_begin()-r.cartesian_end()<<" "<<
                //         std::endl;
                for( ; cit < ait->end() ; ++cit, ++rit) {
                    sum_nom -= *cit * (*rit);
                }
                lamdas[ait-_A.begin()] = sum_nom;
                if(sum_denom==K(0)) {
                    //std::cout<<"div0"<<std::endl;
                    ;
                }
                else {
                    lamda = sum_nom*(1/sum_denom);

                    if(min_plus_not_set && lamda>0) {
                        min_plus=lamda;
                        min_plus_not_set=false;
                    }
                    if(max_minus_not_set && lamda<0) {
                        max_minus=lamda;
                        max_minus_not_set=false;
                    }
                    if(lamda<min_plus && lamda>0) min_plus=lamda;
                    if(lamda>max_minus && lamda<0) max_minus=lamda;
                    //TEST
                    //if(min_plus_not_set && lamda>0){min_plus=lamda;min_plus_not_set=false;mini=ait-_A.begin();}
                    //if(max_minus_not_set && lamda<0){max_minus=lamda;max_minus_not_set=false;
                    //	 maxi=ait-_A.begin();}
                    //if(lamda<min_plus && lamda>0) {min_plus=lamda;mini=ait-_A.begin();}
                    //if(lamda>max_minus && lamda<0) {max_minus=lamda;maxi=ait-_A.begin();}
                }
            }
        } else {//only a few opers no innerprod
            for(typename stdMatrix::iterator ait=_A.begin(); ait<_A.end(); ++ait) {
                typename stdCoeffs::iterator cit;
                cit=ait->begin();
                ++cit;

                NT c_rand_coord = *(cit+rand_coord);
                NT c_rand_coord_prev = *(cit+rand_coord_prev);

                *lamdait = *lamdait
                           + c_rand_coord_prev * (r_prev[rand_coord_prev] - r[rand_coord_prev]);

                if(c_rand_coord==K(0)) {
                    //std::cout<<"div0"<<std::endl;
                    ;
                } else {
                    lamda = (*lamdait) / c_rand_coord;

                    if(min_plus_not_set && lamda>0) {
                        min_plus=lamda;
                        min_plus_not_set=false;
                    }
                    if(max_minus_not_set && lamda<0) {
                        max_minus=lamda;
                        max_minus_not_set=false;
                    }
                    if(lamda<min_plus && lamda>0) min_plus=lamda;
                    if(lamda>max_minus && lamda<0) max_minus=lamda;
                    //TEST
                    //if(min_plus_not_set && lamda>0){min_plus=lamda;min_plus_not_set=false;mini=ait-_A.begin();}
                    //if(max_minus_not_set && lamda<0){max_minus=lamda;max_minus_not_set=false;
                    //	 maxi=ait-_A.begin();}
                    //if(lamda<min_plus && lamda>0) {min_plus=lamda;mini=ait-_A.begin();}
                    //if(lamda>max_minus && lamda<0) {max_minus=lamda;maxi=ait-_A.begin();}
                }
                ++lamdait;
            }
        }
        //std::cout<<"Oresult: "<<mini<<" "<<maxi<<std::endl;
        return std::pair<NT,NT> (min_plus,max_minus);
    }

    //void rotate(){
    //std::cout<<_A<<std::endl;
    //  exit(1);
    //}

private:
    int            _d; //dimension
    stdMatrix      _A; //inequalities
	std::vector<Point>  _sites;
	std::shared_ptr<falconn::LSHNearestNeighborTable<falconn::DenseVector<K>>> hptable;
	std::vector<falconn::DenseVector<K>> data;
	ANNkd_tree* kdTree;
	ANNkd_tree* bdTree;
	Eigen::MatrixXd proj_matrix;
	double _maxDistToBoundary;
	int _k;
	int _l;
    //EXPERIMENTAL
    //Flann_trees    flann_trees; //the (functional) duals of A lifted to answer NN queries
    //defined for every d coordinate
};


// define different kind of polytopes
typedef std::vector<Hyperplane> H_polytope;
typedef H_polytope Polytope;
typedef std::vector<Point> V_polytope;

typedef std::pair<V_polytope,V_polytope> MinkSumPolytope;
typedef std::pair<MinkSumPolytope,bool> MinkSumPolytopeDual;

/* Construct a n-CUBE H-REPRESENTATION*/
Polytope cube(int n, NT lw, NT up) {
    Polytope cube;
    std::vector<NT> origin(n,NT(lw));
    for(int i=0; i<n; ++i) {
        std::vector<NT> normal;
        for(int j=0; j<n; ++j) {
            if(i==j)
                normal.push_back(NT(1));
            else normal.push_back(NT(0));
        }
        Hyperplane h(Point(n,origin.begin(),origin.end()),
                     Direction(n,normal.begin(),normal.end()));
        cube.push_back(h);
    }
    std::vector<NT> apex(n,NT(up));
    for(int i=0; i<n; ++i) {
        //std::cout<<apex[i]<<" ";
        std::vector<NT> normal;
        for(int j=0; j<n; ++j) {
            if(i==j)
                normal.push_back(NT(-1));
            else normal.push_back(NT(0));
        }
        Hyperplane h(Point(n,apex.begin(),apex.end()),
                     Direction(n,normal.begin(),normal.end()));
        cube.push_back(h);
    }
    return cube;
}

/* Construct a n-CUBE V-REPRESENTATION*/
V_polytope Vcube(int n, NT lw, NT up) {
    V_polytope cube;
    for(int k=-1; k<2; k+=2) {
        for(int i=0; i<std::pow(2,n-1); ++i) {
            //bool bytes[sizeof i];
            //std::copy(static_cast<const bool*>(static_cast<const void*>(&i)),
            //          static_cast<const bool*>(static_cast<const void*>(&i)) + sizeof i,
            //          bytes);
            //for(int j=0; j<(sizeof i); ++j)
            boost::dynamic_bitset<> b( n, i );
            std::vector<NT> normal;
            normal.push_back(NT(-1*up*k));
            for (boost::dynamic_bitset<>::size_type j = 0; j < b.size(); ++j) {
                if(b[j]) normal.push_back(NT(1*up));
                else normal.push_back(NT(-1*up));
            }
            //Vector normal_v(n,normal.begin(),normal.end());
            //std::cout<<Vector(n,normal.begin(),normal.end())<<std::endl;
            //std::cout<<Point(n,normal.begin(),normal.end())<<std::endl;
            cube.push_back(Point(n,normal.begin(),normal.end()));
        }
    }
    return cube;
}

/* Construct a n-CROSSPOLYTOPE */
Polytope cross(int n, NT lw, NT up) {
    Polytope cross;
    for(int k=-1; k<2; k+=2) {
        std::vector<NT> vertex;
        vertex.push_back(NT(k*up));
        for(int j=1; j<n; ++j)
            vertex.push_back(NT(0));
        //std::cout<<Point(n,vertex.begin(),vertex.end())<<std::endl;

        for(int i=0; i<std::pow(2,n-1); ++i) {
            //bool bytes[sizeof i];
            //std::copy(static_cast<const bool*>(static_cast<const void*>(&i)),
            //          static_cast<const bool*>(static_cast<const void*>(&i)) + sizeof i,
            //          bytes);
            //for(int j=0; j<(sizeof i); ++j)
            boost::dynamic_bitset<> b( n, i );
            std::vector<NT> normal;
            normal.push_back(NT(-1*k*up));
            for (boost::dynamic_bitset<>::size_type j = 0; j < b.size(); ++j) {
                if(b[j]) normal.push_back(NT(1*up));
                else normal.push_back(NT(-1*up));
            }
            //Vector normal_v(n,normal.begin(),normal.end());
            //std::cout<<Vector(n,normal.begin(),normal.end())<<std::endl;
            Hyperplane h(Point(n,vertex.begin(),vertex.end()),
                         Direction(n,normal.begin(),normal.end()));
            cross.push_back(h);
        }
        //std::cout<<"----"<<std::endl;
    }
    return cross;
}

/* Construct a SKINNY n-CROSSPOLYTOPE */
Polytope cross_skinny(int n, NT lw, NT up) {
    Polytope cross;
    NT sf=pow(2,n);//skinny_factor
    for(int k=-1; k<2; k+=2) {
        std::vector<NT> vertex;
        vertex.push_back(NT(k));
        for(int j=1; j<n; ++j)
            vertex.push_back(NT(0));
        //std::cout<<Point(n,vertex.begin(),vertex.end())<<std::endl;

        for(int i=0; i<std::pow(2,n-1); ++i) {
            //bool bytes[sizeof i];
            //std::copy(static_cast<const bool*>(static_cast<const void*>(&i)),
            //          static_cast<const bool*>(static_cast<const void*>(&i)) + sizeof i,
            //          bytes);
            //for(int j=0; j<(sizeof i); ++j)
            boost::dynamic_bitset<> b( n, i );
            std::vector<NT> normal;
            normal.push_back(NT(-1*k*sf));
            for (boost::dynamic_bitset<>::size_type j = 0; j < b.size(); ++j) {
                if(b[j]) normal.push_back(NT(1));
                else normal.push_back(NT(-1));
            }
            //Vector normal_v(n,normal.begin(),normal.end());
            //std::cout<<Vector(n,normal.begin(),normal.end())<<std::endl;
            Hyperplane h(Point(n,vertex.begin(),vertex.end()),
                         Direction(n,normal.begin(),normal.end()));
            cross.push_back(h);
        }
        //std::cout<<"----"<<std::endl;
    }
    return cross;
}

//SKINNY 2
Polytope cross_skinny2(int n, NT lw, NT up) {
    Polytope cross;
    NT sf=pow(2,n);//skinny_factor
    for(int k=-1; k<2; k+=2) {
        std::vector<NT> vertex;
        vertex.push_back(NT(k*sf));
        for(int j=1; j<n; ++j)
            vertex.push_back(NT(0));
        //std::cout<<Point(n,vertex.begin(),vertex.end())<<std::endl;

        for(int i=0; i<std::pow(2,n-1); ++i) {
            //bool bytes[sizeof i];
            //std::copy(static_cast<const bool*>(static_cast<const void*>(&i)),
            //          static_cast<const bool*>(static_cast<const void*>(&i)) + sizeof i,
            //          bytes);
            //for(int j=0; j<(sizeof i); ++j)
            boost::dynamic_bitset<> b( n, i );
            std::vector<NT> normal;
            normal.push_back(NT(-1*k/sf));
            for (boost::dynamic_bitset<>::size_type j = 0; j < b.size(); ++j) {
                if(b[j]) normal.push_back(NT(1));
                else normal.push_back(NT(-1));
            }
            //Vector normal_v(n,normal.begin(),normal.end());
            //std::cout<<Vector(n,normal.begin(),normal.end())<<std::endl;
            Hyperplane h(Point(n,vertex.begin(),vertex.end()),
                         Direction(n,normal.begin(),normal.end()));
            cross.push_back(h);
        }
        //std::cout<<"----"<<std::endl;
    }
    return cross;
}

/* Construct a n-CROSS V-REPRESENTATION*/
V_polytope Vcross(int n, NT lw, NT up) {
    V_polytope cross;
    for(int i=0; i<n; ++i) {
        std::vector<NT> normal;
        for(int j=0; j<n; ++j) {
            if(i==j)
                normal.push_back(NT(1*up));
            else normal.push_back(NT(0));
        }
        cross.push_back(Point(n,normal.begin(),normal.end()));
        //std::cout<<Point(n,normal.begin(),normal.end())<<std::endl;
    }
    for(int i=0; i<n; ++i) {
        //std::cout<<apex[i]<<" ";
        std::vector<NT> normal;
        for(int j=0; j<n; ++j) {
            if(i==j)
                normal.push_back(NT(-1*up));
            else normal.push_back(NT(0));
        }
        //std::cout<<Point(n,normal.begin(),normal.end())<<std::endl;
        cross.push_back(Point(n,normal.begin(),normal.end()));
    }
    return cross;
}


// contruct a n-ball of radius r centered in the origin
/*
Ball ball(int n, const NT r){

  std::vector<Point> P_ball;
  for(int i=0; i<n; ++i){
		std::vector<NT> coords;
		for(int j=0; j<n; ++j){
			if(i==j)
				coords.push_back(r);
			else coords.push_back(NT(0));
		}
		P_ball.push_back(Point(n,coords.begin(),coords.end()));
	}
	std::vector<NT> extra_coords(n,NT(0));
	extra_coords[0]=NT(-1*r);
	P_ball.push_back(Point(n,extra_coords.begin(),extra_coords.end()));
  Ball B(n,P_ball.begin(),P_ball.end());
	return B;
}
*/

//template <typename T> struct Oracle{
//  sep Sep_Oracle(T &P, Point v);
//};

#endif //POLYTOPES_H
