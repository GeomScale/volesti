// VolEsti (volume computation and sampling library)

// Copyright (c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// VolEsti is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// VolEsti is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with HeaDDaCHe,
// see <http://www.gnu.org/licenses/>.


#ifndef SOLVE_LP_H
#define SOLVE_LP_H

#include <stdio.h>
#include "lp_lib.h"


// compute the chebychev ball of an H-polytope described by a dxd matrix A and  d-dimensional vector b, s.t.: Ax<=b
std::pair<Point,NT> solveLP(Eigen::MatrixXd &A, Eigen::VectorXd &b, int d){

	lprec *lp;
	int Ncol=d+1, *colno = NULL, j, ret = 0, m=A.rows();
	REAL *row = NULL;
	
	lp = make_lp(m, Ncol);
	
	if(lp == NULL)
		ret = 1; /* couldn't construct a new model... */
    
    REAL infinite = get_infinite(lp); /* will return 1.0e30 */
    
    if(ret == 0) {
		/* create space large enough for one row */
		colno = (int *) malloc(Ncol * sizeof(*colno));
		row = (REAL *) malloc(Ncol * sizeof(*row));
		if((colno == NULL) || (row == NULL))
			ret = 2;
	}
	
	
	if(ret == 0) {
		set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */
	}
	int i=0;
	NT sum;
	while(ret==0 & i<m){
		/* construct all rows */
		sum=NT(0);
		for(j=0; j<d; j++){
			colno[j] = j+1;
			row[j] = A(i,j);
			sum+=A(i,j)*A(i,j);
		}
		colno[d] = d+1; /* last column */
		row[d] = std::sqrt(sum);

        /* add the row to lpsolve */
		if(!add_constraintex(lp, d+1, row, colno, LE, b(i))){
			ret = 3;
		}
		i++;
	}
	
	if(ret == 0) {
		set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */


		for (j = 0; j < d; j++) {
            colno[j] = j + 1;
            row[j] = 0;
            set_bounds(lp, j + 1, -infinite, infinite);
        }
        colno[d] = d + 1; /* last column */
        row[d] = 1.0;
        set_bounds(lp, d + 1, 0.0, infinite);

		// set the objective function
		if (!set_obj_fnex(lp, d + 1, row, colno)) {
			ret = 4;
		}
	}
	
	
	if(ret == 0) {
		/* set the object direction to maximize */
		set_maxim(lp);

		/* I only want to see important messages on screen while solving */
		set_verbose(lp, NEUTRAL);

		/* Now let lpsolve calculate a solution */
		ret = solve(lp);
		if(ret == OPTIMAL)
			ret = 0;
		else
			ret = 5;
	}

    std::pair<Point,NT> res;

	if(ret == 0) {
        std::vector<NT> temp_p(d,0);
        get_variables(lp, row);
        for(j = 0; j < d; j++){
			temp_p[j]=NT(row[j]);
	    }
	
	    Point xc( d , temp_p.begin() , temp_p.end() );
	    NT r=NT(get_objective(lp));
	    res = std::pair<Point,NT> (xc,r);
	} else {
        delete_lp(lp);
        std::cout<<"Linear program for the computation of the chebychev ball failed"<<std::endl;
        exit(-1);
    }

    delete_lp(lp);
	return res;
}


// return true if q belongs to the convex hull of the V-polytope described by matrix V
// otherwise return false
bool memLP_Vpoly(Eigen::MatrixXd V, Point q){

	int d=q.dimension();
	lprec *lp;
	int Ncol=d+1, *colno = NULL, j, ret = 0, m=V.rows();
	m++;
	REAL *row = NULL;

	lp = make_lp(m, Ncol);

	if(lp == NULL)
		ret = 1; /* couldn't construct a new model... */

	REAL infinite = get_infinite(lp); /* will return 1.0e30 */

	if(ret == 0) {
		/* create space large enough for one row */
		colno = (int *) malloc(Ncol * sizeof(*colno));
		row = (REAL *) malloc(Ncol * sizeof(*row));
		if((colno == NULL) || (row == NULL))
			ret = 2;
	}

	if(ret == 0) {
		set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */
	}
	int i=0;
	while(ret==0 & i<m-1){
		/* construct all rows */
		for(j=0; j<d; j++){
			colno[j] = j+1;
			row[j] = V(i,j);
		}
		colno[d] = d+1;
		row[d] = -1.0;
		//set_bounds(lp, d, 0.0, infinite);


		/* add the row to lpsolve */
		if(!add_constraintex(lp, d+1, row, colno, LE, 0.0)){
			ret = 3;
		}
		i++;
	}
	if (ret==0){
		for(j=0; j<d; j++){
			colno[j] = j+1; /* last column */
			row[j] = q[j];
		}
		colno[d] = d+1; /* last column */
		row[d] = -1.0;

		/* add the row to lpsolve */
		if(!add_constraintex(lp, d+1, row, colno, LE, 1.0)){
			ret = 3;
		}
	}

	//set the bounds
	if(ret == 0) {
		set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */

        // set the bounds
		for(j=0; j<d; j++){
			colno[j] = j+1; /* j_th column */
			row[j] = q[j];
			set_bounds(lp, j+1, -infinite, infinite);
		}
		colno[d] = d+1; /* last column */
		row[d] = -1.0;
		set_bounds(lp, d+1, -infinite, infinite);

        // set the objective function
		if(!set_obj_fnex(lp, d+1, row, colno)){
			ret = 4;
		}

	}

	if(ret == 0) {
		/* set the object direction to maximize */
		set_maxim(lp);

		/* I only want to see important messages on screen while solving */
		set_verbose(lp, NEUTRAL);

		/* Now let lpsolve calculate a solution */
		ret = solve(lp);
		if(ret == OPTIMAL)
			ret = 0;
		else
			ret = 5;
	}

	if (ret == 0) {
        NT r = NT(get_objective(lp));
        delete_lp(lp);
        if(r>0.0){
            return false;
        }
        return true;
	}

    std::cout<<"Linear Program for the membership failed"<<std::endl;
    exit(-1);
    return 0;
}


// compute the intersection of a ray with a V-polytope
// if maxi is true compute positive lambda, when the ray is p + lambda \codt v
// otherwise compute the negative lambda
NT intersect_line_Vpoly(Eigen::MatrixXd V, Point &p, Point &v, bool maxi){

	int d=v.dimension();
	lprec *lp;
	int m=V.rows();
	m++;
	int Ncol=m, *colno = NULL, j, ret = 0;
	REAL *row = NULL;
	NT res;

	lp = make_lp(d+1, Ncol);

	if(lp == NULL)
		ret = 1; /* couldn't construct a new model... */

	REAL infinite = get_infinite(lp); /* will return 1.0e30 */

	if(ret == 0) {
		/* create space large enough for one row */
		colno = (int *) malloc(Ncol * sizeof(*colno));
		row = (REAL *) malloc(Ncol * sizeof(*row));
		if((colno == NULL) || (row == NULL))
			ret = 2;
	}

	if(ret == 0) {
		set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */
	}
	int i=0;

	while(ret==0 & i<d){
		/* construct all rows  */
		for(j=0; j<m-1; j++){
			colno[j] = j+1; /* j_th column */
			row[j] = V(j,i);
		}
		colno[m-1] = m; /* last column */
		row[m-1] = v[i];

		/* add the row to lpsolve */
		if(!add_constraintex(lp, m, row, colno, EQ, p[i])){
			ret = 3;
		}
		i++;
	}

	if(ret==0){
		for(j=0; j<m-1; j++){
			colno[j] = j+1; /* j_th column */
			row[j] = 1.0;
		}
		colno[m-1] = m; /* last column */
		row[m-1] = 0.0;

		/* add the row to lpsolve */
		if(!add_constraintex(lp, m, row, colno, EQ, 1.0)){
			ret = 3;
		}
	}

	//set the bounds
	if(ret == 0) {
		set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */

		// set the objective function
		for(j=0; j<m-1; j++){
			colno[j] = j+1; /* j_th column */
			row[j] = 0;
		}
		colno[m - 1] =m; /* last column */
		row[m-1] = 1.0;
		set_bounds(lp, m, -infinite, infinite);

		// set objective function
		if(!set_obj_fnex(lp, m, row, colno)){
			ret = 4;
		}

	}

	if(ret == 0) {

		if(maxi) {  /* set the object direction to maximize */
			set_maxim(lp);
		}else{      /* set the object direction to minimize */
			set_minim(lp);
		}
		set_verbose(lp, NEUTRAL);

		/* Now let lpsolve calculate a solution */
		ret = solve(lp);
		if(ret == OPTIMAL)
			ret = 0;
		else
			ret = 5;
	}

	if(ret==0) {
		res = NT(-get_objective(lp));
		delete_lp(lp);
	    return res;

	}
	std::cout<<"Linear Program for hit and run failed"<<std::endl;
	exit(-1);
	return 0;
}

#endif
