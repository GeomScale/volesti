// Copyright (c) 2010 ETH Zurich (Switzerland)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
// 
// $URL: $
// $Id:  $
// 
//
// Author(s)     : Christian Helbling

#ifndef CGAL_EXTREME_POINTS_D_H
#define CGAL_EXTREME_POINTS_D_H
#include <CGAL/config.h>
#include <CGAL/basic.h>
#include <vector>
#include <algorithm>
#include <iterator>

#include <CGAL/QP_functions.h>
#include <CGAL/QP_solution.h>
#include <CGAL/Quotient.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Coercion_traits.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Extreme_points_traits_d.h>
#include <CGAL/Extreme_points_options_d.h>

#include <boost/functional.hpp>
#include <boost/function.hpp>
#include <boost/iterator.hpp>

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpzf.h>
typedef CGAL::Gmpzf ETF;
typedef CGAL::Gmpz  ETI;
#else
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ETF;
typedef CGAL::MP_Float ETI;
#endif

namespace CGAL {

// templates for finding an appropriate exact type..
namespace internal {
    
    template <class T, class Is_exact_tag, class Algebraic_category_tag>
    struct Exact_type_base {
        typedef ETF Type; // exact floating type as default
    };
    
    template <class T, class Algebraic_category_tag>
    struct Exact_type_base<T, CGAL::Tag_true, Algebraic_category_tag> {
        typedef T Type; // T is already an exact type
    };
    
    template <class T>
    struct Exact_type_base<T, CGAL::Tag_false, CGAL::Field_tag> {
        typedef ETI Type; // T is some integer type..
    };
    
    template <class T>
    struct Exact_type {
        typedef Algebraic_structure_traits<T> AST;
        typedef typename Exact_type_base<T, typename AST::Is_exact,
            typename AST::Algebraic_category>::Type ET;
//         make sure the chosen exact type is also compliant with T
        typedef typename Coercion_traits<T,ET>::Type Type;
    };
}

template <class InputIterator, class OutputIterator>
OutputIterator
extreme_points_d_dula_helgason(InputIterator first, InputIterator beyond,
                    OutputIterator  result);

template <class InputIterator, class OutputIterator>
OutputIterator
extreme_points_d_simple(InputIterator first, InputIterator beyond,
                                OutputIterator  result);

enum Extreme_point_classification {INTERNAL_POINT=-1, EXTREME_POINT=0,
                                   EXTERNAL_POINT=1 };

template <class Traits>
class Extreme_points_d {
    public:
        // types
        typedef typename Traits::Point                          Point;
        typedef typename Traits::Less_lexicographically         Less_lexicographically;
        typedef typename Traits::RT                             RT;

    private:
        // the exact type which is used in the QP-Solver
        typedef typename CGAL::internal::Exact_type<RT>::Type   ET;
    
    protected:
        int dim;
        
        Extreme_points_options_d ep_options_;
        
        // points inserted after the last extreme point computation
        std::vector<Point> new_points;
        
        // result of the last extreme point computation
        std::vector<Point> extreme_points;
        
        // update the extreme_points vector if we have new points
        void update();
        
    public:
        Extreme_points_d(int d, Extreme_points_options_d ep_options = 
                                Extreme_points_options_d())
            : dim(d), ep_options_(ep_options) {}
        
        int dimension() {
            return dim;
        }
        
        void clear() {
            new_points.clear();
            extreme_points.clear();
        }
        
        void insert(const Point x) {
            CGAL_precondition_msg(x.dimension() == dim,
                                  "Invalid dimension of inserted point.");
            new_points.push_back(x);
        }
        
        template <typename InputIterator>
        void insert(InputIterator first, InputIterator beyond) {
            while (first != beyond) insert(*first++);
        }
        
        // generates the extreme points of the given point set
        template <class OutputIterator>
        OutputIterator
        get_extreme_points(OutputIterator  result) {
            update();
            return std::copy(extreme_points.begin(), extreme_points.end(),
                             result);
        }
        
        // classifies the given point as internal, extreme or external point
        // if is_input_point is set, p is assumed to be either an
        // extreme point or an internal point
        enum Extreme_point_classification classify(Point p,
                                                   bool is_input_point=false);
};

namespace internal {
    // calculates the inner product of a d-dimensional point
    // (using homogeneous coordinates)
    // with a d-vector (Cartesian coordinates of type ET)
    template <class Point, class ET>
    Quotient<ET>
    inner_product_pv(Point p, std::vector<ET> &u) {
        ET r=0;
        for (int i=0;i<p.dimension();++i)
            r+=u[i]*ET(p.homogeneous(i));
        return Quotient<ET> (r,p.homogeneous(p.dimension()));
    }
    
    // function to solve the LP that tests whether a point is in the
    // convex hull of other points; the type ET is an exact type used
    // for the internal computations
    template <class Point, class ET, class RandomAccessIterator, class Traits>
    CGAL::Quadratic_program_solution<ET>
    solve_convex_hull_containment_lp (const Point p,
                                      RandomAccessIterator begin,
                                      RandomAccessIterator end,
                                      const ET &et_dummy,
                                      const Traits &ep_traits,
                                      const Quadratic_program_options &qp_options) {
        typedef typename Traits::Homogeneous_begin  Homogeneous_begin;
        
        // Constraint matrix type: A[j][i] is the i-th homogeneous coordinate
        // of a_j
        typedef boost::transform_iterator <Homogeneous_begin,
                                           RandomAccessIterator> A_it;
        
        // Right-hand side type: b[i] is the i-th homogeneous coordinate of p
        typedef typename Homogeneous_begin::result_type B_it;
        
        // Relation type ("=")
        typedef CGAL::Const_oneset_iterator<CGAL::Comparison_result> R_it;
        
        // input number type
        typedef typename Traits::RT RT;
        
        // Linear objective function type (c=0: we only test feasibility)
        typedef CGAL::Const_oneset_iterator<RT> C_it;
        
        // the nonnegative linear program type
        typedef
            CGAL::Nonnegative_linear_program_from_iterators<A_it, B_it, R_it,
                                                            C_it> LP;
        
        LP lp (end-begin,              // number of variables
               p.dimension()+1,        // number of constraints
               A_it (begin),           // A
               Homogeneous_begin()(p), // b
               R_it (CGAL::EQUAL),     // r
               C_it (0));              // c
        
        return CGAL::solve_nonnegative_linear_program (lp, ET(0), qp_options);
    }
    
    template <class Point, class RandomAccessIterator, class ET, class Traits>
    bool is_in_convex_hull (const Point p,
                            RandomAccessIterator begin,
                            RandomAccessIterator end,
                            const ET& et_dummy,
                            const Traits &ep_traits,
                            const Quadratic_program_options &qp_options) {
        CGAL::Quadratic_program_solution<ET> s =
            solve_convex_hull_containment_lp (p, begin, end, et_dummy,
                                              ep_traits, qp_options);
        return !s.is_infeasible();
    }
    
} // namespace internal

template <class Traits>
void Extreme_points_d<Traits>::update() {
    // only do something if we have new points
    if (!new_points.empty()) {
        // chose which frame-algorithm to run
        Extreme_point_algorithm_d algo = ep_options_.get_algorithm();
        if (algo==EP_CHOOSE_APPROPRIATE) {
            // our simple heuristic is just the ratio of new points to the old
            // extreme points
            // if we only have few new points, then we expect most of the 
            // points we run the algorithm on to stay extreme, so we have
            // a big extreme point ratio in which case the simple algorithm 
            // performs slightly better
            
            if (4 * new_points.size() < extreme_points.size())
                algo=EP_SIMPLE;
            else
                algo=EP_DULA_HELGASON;
        }
        
        
        // temporarily add extreme_points to new_points
        new_points.insert(new_points.begin(), extreme_points.begin(),
                          extreme_points.end());
        
        // delete extreme_points as we fill it with the new extreme points
        extreme_points.clear();

        
        // run the frame algorithm
        switch (algo) {
            case EP_SIMPLE:
                extreme_points_d_dula_helgason(new_points.begin(), new_points.end(),
                                               std::back_inserter(extreme_points),
                                               Traits(),
                                               ep_options_.get_qp_options() );
                break;
            case EP_DULA_HELGASON:
            default:
                extreme_points_d_dula_helgason(new_points.begin(), new_points.end(),
                                               std::back_inserter(extreme_points),
                                               Traits(),
                                               ep_options_.get_qp_options() );
                break;
                
        };
        
        // sort extreme points as we want to use this array for binary search in the classify function
        sort(extreme_points.begin(),extreme_points.end(),
             Less_lexicographically());
        
        // clear the new_points array
        new_points.clear();
    }
}

template <class Traits>
enum Extreme_point_classification
Extreme_points_d<Traits>::classify(Point p, bool is_input_point) {
    update();
    
    if  (std::binary_search(extreme_points.begin(),extreme_points.end(),p,
                            Less_lexicographically())) {
        // p is an extreme point
        return EXTREME_POINT;
    } else if (is_input_point || 
               CGAL::internal::is_in_convex_hull(p,
                                                 extreme_points.begin(),
                                                 extreme_points.end(),
                                                 ET(0),
                                                 Traits(),
                                                 ep_options_.get_qp_options())) {
        // p is an internal point
        return INTERNAL_POINT;
    } else {
        // p is some external point
        return EXTERNAL_POINT;
    }
}

/**
Algorithm by Dulá and Helgason from "A new procedure for identifying the frame
of the convex hull of a finite collection of points in multidimensional space",
1996 
as described as OS_Frame01 in "Competing Output-Sensitive Frame Algorithms"
(J.H. Dulá, F.J. López)

Own additions:
as we don't assume the input data to be in general position:
* duplicated points are handled as one point
* when chosing the new extreme point in the case that the tested point is
not in the convex hull of the frame elements found so far we consider
lexicographical ordering in case multiple points have the same scalar
product
**/
template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
extreme_points_d_dula_helgason(InputIterator first, InputIterator beyond,
                    OutputIterator  result, const Traits &ep_traits,
                    const Quadratic_program_options &qp_options =
                        Quadratic_program_options() ) {
    typedef typename Traits::Point                      Point;
    typedef typename Traits::Less_lexicographically     Less_lexicographically;
    typedef typename Traits::RT                         RT;

    // the exact type which is used in the QP-Solver
    typedef typename CGAL::internal::Exact_type<RT>::Type
                                                        ET;
    typedef CGAL::Quadratic_program_solution<ET>        QP_Solution;
    
    std::vector<Point> points(first,beyond); // all input points
    
    int n=points.size();
    if (n==0)
        return result; // no points -> nothing to do..
    
    int d=points[0].dimension();
    
    // move the lexicographically smallest point to the front
    std::nth_element(points.begin(), points.begin(), points.end(),
                     Less_lexicographically());
    
    // now shuffle the points for better performance
    // let the first point stay, as we begin with that one..
    std::random_shuffle(points.begin()+1,points.end());
    
    std::vector<bool> candidate(n, true);
    
    // initialize f with the first point
    // which has the lexicographically smallest coordinates
    CGAL_assertion(n>0);
    std::vector<Point> f;
    f.push_back(points[0]);
    *result++=points[0]; // already the first extreme point
    candidate[0]=false;
    
    // start at point 1 as we already took point 0
    for (int j=1;j<n;++j) {
        if (!candidate[j])
            continue; // point[j] already in conv(f)
        while (1) {
            QP_Solution s =
                CGAL::internal::solve_convex_hull_containment_lp(
                    points[j], f.begin(), f.end(), ET(), ep_traits, qp_options);
            
            if (s.is_infeasible()) {
                // points[j] \notin conv(f)
                
                // get certificate (which is the hyperplane pi which separates
                // f from points[j])
                typename QP_Solution::Infeasibility_certificate_iterator ifc =
                    s.infeasibility_certificate_begin();
                
                // store the separating hyperplane in a vector<ET> called sh
                std::vector<ET> sh(d+1);
                for (int i=0;i<=d;++i)
                    sh[i]=ifc[i];
                
                // we are looking for the point farthest away from the
                // separating hyperplane sh on the side of points[j]
                Quotient<ET> m = CGAL::internal::inner_product_pv(points[j],
                                                                  sh);
                
                int k=j;
                
                for (int i=j+1;i<n;++i) {
                    if (!candidate[i])
                        continue;
                    
                    Quotient<ET> x=CGAL::internal::inner_product_pv(points[i],
                                                                    sh);
                    
                    if (x<m ||
                        x==m && Less_lexicographically()(points[i],points[k]))
                    {
                        m=x;
                        k=i;
                    }
                }
                
                *result++=points[k];
                candidate[k]=false;
                f.push_back(points[k]);
                if (k==j) // now points[j] \in conv(f) for sure
                    break;
            } else {
                // points[j] \in conv(f)
                candidate[j]=false; // not strictly nessessary
                break;
            }
        }
    }
    return result;
}

template <class InputIterator, class OutputIterator>
OutputIterator
extreme_points_d_dula_helgason(InputIterator first, InputIterator beyond,
                               OutputIterator  result) {
    typedef std::iterator_traits<InputIterator> ITraits;
    typedef typename ITraits::value_type        Point;
    typedef Extreme_points_traits_d<Point>      Traits;
    return extreme_points_d_dula_helgason(first, beyond, result, Traits());
}

// generates the extreme points of the given point set
// simple (naive) algorithm
// no general position assumptions
template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
extreme_points_d_simple(InputIterator first, InputIterator beyond,
                        OutputIterator  result, const Traits &ep_traits,
                        const Quadratic_program_options &qp_options =
                            Quadratic_program_options() ) {
    typedef typename Traits::Point                      Point;
    typedef typename Traits::Less_lexicographically     Less_lexicographically;
    typedef typename Traits::RT                         RT;
    
    // the exact type which is used in the QP-Solver
    typedef typename CGAL::internal::Exact_type<RT>::Type
                                                        ET;
    typedef CGAL::Quadratic_program_solution<ET>        QP_Solution;
    
    std::vector<Point> points_in(first,beyond); // all input points
    
    int n=points_in.size();
    
    if (n==0)
        return result; // nothing to do
    
    // remove duplicate points
    sort(points_in.begin(), points_in.end(), Less_lexicographically());
    
    std::vector<Point> points; // pointers to all distinct points
    points.push_back(points_in[0]);
    for (int i=1;i<n;++i) {
        // don't take the duplicated ones more than once
        if (points_in[i-1] != points_in[i]) {
            points.push_back(points_in[i]);
        }
    }
    // adjust n
    n=points.size();

    // shuffle again so that the LP-solver runs faster
    std::random_shuffle(points.begin(),points.end());
    
    // test all points for being in the convex hull of the others
    // using LP
    for (int i=0;i<n;++i) {
        // test i-th point
        std::swap(points[i],points[0]); // move point to test at postition 0
        
        if (!CGAL::internal::is_in_convex_hull(points[0],
                                               points.begin()+1,
                                               points.end(),
                                               ET(0),
                                               Traits(),
                                               qp_options)) {
            *result++=points[0];
        }
        std::swap(points[i],points[0]); // move test point back
    }
    return result;
}

template <class InputIterator, class OutputIterator>
OutputIterator
extreme_points_d_simple(InputIterator first, InputIterator beyond,
                        OutputIterator  result) {
    typedef std::iterator_traits<InputIterator> ITraits;
    typedef typename ITraits::value_type        Point;
    typedef Extreme_points_traits_d<Point>      Traits;
    return extreme_points_d_simple(first, beyond, result, Traits());
}

template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
extreme_points_d(InputIterator first, InputIterator beyond,
                 OutputIterator  result, const Traits &ep_traits,
                 const Quadratic_program_options &qp_options =
                     Quadratic_program_options() ) {
    return extreme_points_d_dula_helgason(first, beyond, result, ep_traits, qp_options);
}

template <class InputIterator, class OutputIterator>
OutputIterator
extreme_points_d(InputIterator first, InputIterator beyond,
                 OutputIterator  result) {
    typedef std::iterator_traits<InputIterator> ITraits;
    typedef typename ITraits::value_type        Point;
    typedef Extreme_points_traits_d<Point>      Traits;
    return extreme_points_d(first, beyond, result, Traits());
}

} //namespace CGAL

#endif // CGAL_EXTREME_POINTS_D_H
