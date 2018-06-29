#include <stdio.h>
#include "lp_lib.h"


template <typename K>
std::pair<Point,double> solveLP(std::vector<std::vector<K> > A, int d){
	
	lprec *lp;
	int Ncol=d+1, *colno = NULL, j, ret = 0, m=A.size();
	REAL *row = NULL;
	
	lp = make_lp(m, Ncol);
	
	if(lp == NULL)
		ret = 1; /* couldn't construct a new model... */
    
    REAL infinite = get_infinite(lp); /* will return 1.0e30 */
    
    if(ret == 0) {
    /* let us name our variables. Not required, but can be useful for debugging */
		
		//for (int i=0; i<d; i++){
			//set_col_name(lp, i+1, "xi");
			//set_col_name(lp, 2, "y");
		//}

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
	K sum;
	while(ret==0 & i<m){
		/* construct all rows (120 x + 210 y <= 15000) */
		sum=K(0);
		for(j=0; j<d; j++){
			colno[j] = j+1; /* j_th column */
			row[j] = A[i][j+1];
			sum += A[i][j+1]*A[i][j+1];
		}
		colno[d] = d+1; /* last column */
		row[d] = std::sqrt(sum);
		//set_bounds(lp, d, 0.0, infinite);
		

    /* add the row to lpsolve */
		if(!add_constraintex(lp, d+1, row, colno, LE, A[i][0])){
			ret = 3;
		}
		i++;
	}
	
	if(ret == 0) {
		set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */
		
		// set the objective function
		for(j=0; j<d; j++){
			colno[j] = j+1; /* j_th column */
			row[j] = 0;
			set_bounds(lp, j+1, -infinite, infinite);
		}
		colno[d] = d+1; /* last column */
		row[d] = 1.0;
		set_bounds(lp, d+1, 0.0, infinite);
		if(!set_obj_fnex(lp, d+1, row, colno)){
			ret = 4;
		}
		
	}
	
	
	if(ret == 0) {
		/* set the object direction to maximize */
		set_maxim(lp);

		/* just out of curioucity, now show the model in lp format on screen */
		/* this only works if this is a console application. If not, use write_lp and a filename */
		//write_LP(lp, stdout);
		/* write_lp(lp, "model.lp"); */

		/* I only want to see important messages on screen while solving */
		set_verbose(lp, NEUTRAL);

		/* Now let lpsolve calculate a solution */
		ret = solve(lp);
		if(ret == OPTIMAL)
			ret = 0;
		else
			ret = 5;
	}
	
	//if(ret == 0) {
	std::vector<double> temp_p(d,0);
	get_variables(lp, row);
	for(j = 0; j < d; j++){
	//printf("%s: %f\n", get_col_name(lp, j + 1), row[j]);
		temp_p[j]=double(row[j]);
	}
	
	Point xc( d , temp_p.begin() , temp_p.end() );
	double r=double(get_objective(lp));
	
			
	return std::pair<Point,double> (xc,r);

	
}
