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


#ifndef PERMUTAEDRON_ORACLES_H
#define PERMUTAEDRON_ORACLES_H

template <class MT, class VT, class Point, typename NT>
NT intersect_line_permutaedron(MT &T, MT &A, VT &b, MT &Aeq, MT &beq, Point &r, Point &v,  NT *conv_comb, NT *row, int *colno) {

    int d = v.dimension(), m = A.rows(), i, Ncol = T.cols()+1, j, Nrows, k = T.cols();
    lprec *lp;
    NT res;
    Nrows = m + d;

    try
    {
        lp = make_lp(Nrows, Ncol);
        if(lp == NULL) throw false;
    }
    catch (bool e) {
#ifdef VOLESTI_DEBUG
        std::cout<<"Could not construct Linear Program for ray-shooting "<<e<<std::endl;
#endif
        return -1.0;
    }

    REAL infinite = get_infinite(lp); /* will return 1.0e30 */

    set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */

    for (i=0; i<d; i++){
        /* construct all rows  */
        for(j=0; j<k; j++){
            colno[j] = j+1; /* j_th column */
            row[j] = T(i,j);
        }
        colno[k] = k + 1; /* last column */
        row[k] = v[i];

        /* add the row to lpsolve */
        try {
            if(!add_constraintex(lp, Ncol, row, colno, EQ, r[i])) throw false;
        }
        catch (bool e)
        {
#ifdef VOLESTI_DEBUG
            std::cout<<"Could not construct constaints for the Linear Program for ray-shooting "<<e<<std::endl;
#endif
            return -1.0;
        }

    }

    for (i=0; i<m; i++){
        /* construct all rows  */
        for(j=0; j<k; j++){
            //colno[j] = j+1; /* j_th column */
            row[j] = A(i,j);
        }
        //colno[k] = k + 1; /* last column */
        row[k] = 0.0;

        /* add the row to lpsolve */
        try {
            if(!add_constraintex(lp, Ncol, row, colno, EQ, b(i))) throw false;
        }
        catch (bool e)
        {
#ifdef VOLESTI_DEBUG
            std::cout<<"Could not construct constaints for the Linear Program for ray-shooting "<<e<<std::endl;
#endif
            return -1.0;
        }

    }

    //set the bounds
    set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */

    // set the objective function
    for(j=0; j<Ncol-1; j++){
        set_bounds(lp, j + 1, 0.0, infinite);
        //colno[j] = j+1; /* j_th column */
        row[j] = 0;
    }
    set_bounds(lp, Ncol, -infinite, infinite);
    row[Ncol - 1] = 1.0;

    // set objective function
    try
    {
        if(!set_obj_fnex(lp, Ncol, row, colno)) throw false;
    }
    catch (bool e)
    {
#ifdef VOLESTI_DEBUG
        std::cout<<"Could not construct objective function for the Linear Program for ray-shooting "<<e<<std::endl;
#endif
        return -1.0;
    }

    /* set the object direction to minimize*/
    set_minim(lp);
    //set_verbose(lp, NEUTRAL);

    /* Now let lpsolve calculate a solution */
    try
    {
        if (solve(lp) != OPTIMAL) throw false;
    }
    catch (bool e)
    {
#ifdef VOLESTI_DEBUG
        std::cout<<"Could not solve the Linear Program for ray-shooting "<<e<<std::endl;
#endif
        return -1.0;
    }

    res = NT(-get_objective(lp));
    //std::cout<<"res = "<<res<<std::endl;
    get_variables(lp, conv_comb);
    delete_lp(lp);
    return res;

}


template <class MT, class VT, class Point, typename NT>
std::pair<NT,NT> intersect_double_line_permutaedron2(MT &T, MT &A, VT &b, Point &r, Point &v,  NT *conv_comb, NT *row, int *colno) {

    int d = v.dimension(), m = A.rows(), i, Ncol = T.cols()+1, j, Nrows, k = T.cols();
    lprec *lp;
    NT res;
    Nrows = m + d;
    std::pair<NT,NT> exc_res(0,0);

    try
    {
        lp = make_lp(Nrows, Ncol);
        if(lp == NULL) throw false;
    }
    catch (bool e) {
#ifdef VOLESTI_DEBUG
        std::cout<<"Could not construct Linear Program for ray-shooting "<<e<<std::endl;
#endif
        return exc_res;
    }

    REAL infinite = get_infinite(lp); /* will return 1.0e30 */

    set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */

    for (i=0; i<d; i++){
        /* construct all rows  */
        for(j=0; j<k; j++){
            colno[j] = j+1; /* j_th column */
            row[j] = T(i,j);
        }
        colno[k] = k + 1; /* last column */
        row[k] = v[i];

        /* add the row to lpsolve */
        try {
            if(!add_constraintex(lp, Ncol, row, colno, EQ, r[i])) throw false;
        }
        catch (bool e)
        {
#ifdef VOLESTI_DEBUG
            std::cout<<"Could not construct constaints for the Linear Program for ray-shooting "<<e<<std::endl;
#endif
            return exc_res;
        }

    }

    for (i=0; i<m; i++){
        /* construct all rows  */
        for(j=0; j<k; j++){
            //colno[j] = j+1; /* j_th column */
            row[j] = A(i,j);
        }
        //colno[k] = k + 1; /* last column */
        row[k] = 0.0;

        /* add the row to lpsolve */
        try {
            if(!add_constraintex(lp, Ncol, row, colno, EQ, b(i))) throw false;
        }
        catch (bool e)
        {
#ifdef VOLESTI_DEBUG
            std::cout<<"Could not construct constaints for the Linear Program for ray-shooting "<<e<<std::endl;
#endif
            return exc_res;
        }

    }

    //set the bounds
    set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */

    // set the objective function
    for(j=0; j<Ncol-1; j++){
        set_bounds(lp, j + 1, 0.0, infinite);
        //colno[j] = j+1; /* j_th column */
        row[j] = 0;
    }
    set_bounds(lp, Ncol, -infinite, infinite);
    row[Ncol - 1] = 1.0;

    // set objective function
    try
    {
        if(!set_obj_fnex(lp, Ncol, row, colno)) throw false;
    }
    catch (bool e)
    {
#ifdef VOLESTI_DEBUG
        std::cout<<"Could not construct objective function for the Linear Program for ray-shooting "<<e<<std::endl;
#endif
        return exc_res;
    }

    std::pair<NT,NT> res_pair;
    set_verbose(lp, NEUTRAL);

    set_maxim(lp);
    //solve(lp);
    if (solve(lp) != OPTIMAL){
        //delete_lp(lp);
        std::cout<<"not feasible1"<<std::endl;
    }

    res_pair.second = NT(-get_objective(lp));

    set_minim(lp);
    //solve(lp);
    if (solve(lp) != OPTIMAL){
        //delete_lp(lp);
        std::cout<<"not feasible2"<<std::endl;
    }
    res_pair.first = NT(-get_objective(lp));
    get_variables(lp, conv_comb);
    for (int l = 0; l < Ncol; ++l) {
        //std::cout<<"convcomb["<<l<<"] = "<<conv_comb[l]<<std::endl;
    }
    //NT sum;

    delete_lp(lp);

    //std::cout<<"l1 = "<<res_pair.first<<" l2 = "<<res_pair.second<<std::endl;
    return res_pair;

}


template <class MT, class VT, class Point, typename NT>
std::pair<NT,NT> intersect_double_line_permutaedron(MT &T, MT &A, VT &b, MT &Aeq, VT &beq, Point &r, Point &v,  NT *conv_comb, NT *row, int *colno) {

    int d = v.dimension(), m = Aeq.rows(), i, Ncol = T.cols()+1, j, Nrows, k = T.cols();
    lprec *lp;
    NT res;
    Nrows = m + d + d*d;
    std::pair<NT,NT> exc_res(0,0);

    try
    {
        lp = make_lp(Nrows, Ncol);
        if(lp == NULL) throw false;
    }
    catch (bool e) {
#ifdef VOLESTI_DEBUG
        std::cout<<"Could not construct Linear Program for ray-shooting "<<e<<std::endl;
#endif
        return exc_res;
    }

    REAL infinite = get_infinite(lp); /* will return 1.0e30 */

    set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */

    for (i=0; i<d; i++){
        /* construct all rows  */
        for(j=0; j<k; j++){
            colno[j] = j+1; /* j_th column */
            row[j] = T(i,j);
        }
        colno[k] = k + 1; /* last column */
        row[k] = v[i];

        /* add the row to lpsolve */
        try {
            if(!add_constraintex(lp, Ncol, row, colno, EQ, r[i])) throw false;
        }
        catch (bool e)
        {
#ifdef VOLESTI_DEBUG
            std::cout<<"Could not construct constaints for the Linear Program for ray-shooting "<<e<<std::endl;
#endif
            return exc_res;
        }

    }

    for (i=0; i<m; i++){
        /* construct all rows  */
        for(j=0; j<k; j++){
            //colno[j] = j+1; /* j_th column */
            row[j] = Aeq(i,j);
        }
        //colno[k] = k + 1; /* last column */
        row[k] = 0.0;

        /* add the row to lpsolve */
        try {
            if(!add_constraintex(lp, Ncol, row, colno, EQ, beq(i))) throw false;
        }
        catch (bool e)
        {
#ifdef VOLESTI_DEBUG
            std::cout<<"Could not construct constaints for the Linear Program for ray-shooting "<<e<<std::endl;
#endif
            return exc_res;
        }

    }

    int count =0;
    for (i=0; i<d*d; i++){
        /* construct all rows  */
        for(j=0; j<k; j++){
            row[j] = A(i,j);
        }
        count++;
        //colno[k] = k + 1; /* last column */
        row[k] = 0.0;

        /* add the row to lpsolve */
        try {
            if(!add_constraintex(lp, Ncol, row, colno, LE, b(i))) throw false;
        }
        catch (bool e)
        {
#ifdef VOLESTI_DEBUG
            std::cout<<"Could not construct constaints for the Linear Program for ray-shooting "<<e<<std::endl;
#endif
            return exc_res;
        }

    }

    //set the bounds
    set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */

    // set the objective function
    for(j=0; j<Ncol-1; j++){
        set_bounds(lp, j + 1, -infinite, infinite);
        //colno[j] = j+1; /* j_th column */
        row[j] = 0;
    }
    set_bounds(lp, Ncol, -infinite, infinite);
    row[Ncol - 1] = 1.0;

    // set objective function
    try
    {
        if(!set_obj_fnex(lp, Ncol, row, colno)) throw false;
    }
    catch (bool e)
    {
#ifdef VOLESTI_DEBUG
        std::cout<<"Could not construct objective function for the Linear Program for ray-shooting "<<e<<std::endl;
#endif
        return exc_res;
    }

    std::pair<NT,NT> res_pair;
    set_verbose(lp, NEUTRAL);

    set_maxim(lp);
    //solve(lp);
    if (solve(lp) != OPTIMAL){
        //delete_lp(lp);
        std::cout<<"not feasible1"<<std::endl;
    }

    res_pair.second = NT(-get_objective(lp));

    set_minim(lp);
    //solve(lp);
    if (solve(lp) != OPTIMAL){
        //delete_lp(lp);
        std::cout<<"not feasible2"<<std::endl;
    }
    res_pair.first = NT(-get_objective(lp));
    get_variables(lp, conv_comb);
    for (int l = 0; l < Ncol; ++l) {
        //std::cout<<"convcomb["<<l<<"] = "<<conv_comb[l]<<std::endl;
    }
    //NT sum;

    delete_lp(lp);

    //std::cout<<"l1 = "<<res_pair.first<<" l2 = "<<res_pair.second<<std::endl;
    return res_pair;

}



template <class MT, class VT, class Point, typename NT2>
bool memLP_permutaedron(MT &T, MT &A, VT &b, MT &Aeq, VT &beq, Point &q, NT2 *mem_comb) {

    typedef typename Point::FT NT;
    int d=q.dimension(), m = Aeq.rows(), k = T.cols(), m2 = A.rows();
    lprec *lp;
    int Ncol=k, *colno = NULL, j, i, Nrows = d+m+m2;

    //REAL *row = NULL;

    try
    {
        lp = make_lp(Nrows, Ncol);
        if(lp == NULL) throw false;
    }
    catch (bool e) {
#ifdef VOLESTI_DEBUG
        std::cout<<"Could not construct Linear Program for membership "<<e<<std::endl;
#endif
        return false;
    }

    REAL infinite = get_infinite(lp); /* will return 1.0e30 */

    try
    {
        colno = (int *) malloc(Ncol * sizeof(*colno));
        //row = (REAL *) malloc(Ncol * sizeof(*row));
    }
    catch (std::exception &e)
    {
#ifdef VOLESTI_DEBUG
        std::cout<<"Linear Program for membership failed "<<e.what()<<std::endl;
#endif
        return false;
    }

    set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */

    for (i = 0;  i< d; ++i) {
        /* construct all rows */
        for(j=0; j<k; j++){
            colno[j] = j+1;
            mem_comb[j] = T(i,j);
        }

        /* add the row to lpsolve */
        try {
            if(!add_constraintex(lp, k, mem_comb, colno, EQ, q[i])) throw false;
        }
        catch (bool e)
        {
#ifdef VOLESTI_DEBUG
            std::cout<<"Could not construct constaints for the Linear Program for membership "<<e<<std::endl;
#endif
            return false;
        }
    }

    for (i = 0;  i< m; ++i) {
        /* construct all rows */
        for(j=0; j<k; j++){
            //colno[j] = j+1;
            mem_comb[j] = Aeq(i,j);
        }

        /* add the row to lpsolve */
        try {
            if(!add_constraintex(lp, k, mem_comb, colno, EQ, beq(i))) throw false;
        }
        catch (bool e)
        {
#ifdef VOLESTI_DEBUG
            std::cout<<"Could not construct constaints for the Linear Program for membership "<<e<<std::endl;
#endif
            return false;
        }
    }

    for (i = 0;  i< m2; ++i) {
        /* construct all rows */
        for(j=0; j<k; j++){
            //colno[j] = j+1;
            mem_comb[j] = A(i,j);
        }

        /* add the row to lpsolve */
        try {
            if(!add_constraintex(lp, k, mem_comb, colno, LE, b(i))) throw false;
        }
        catch (bool e)
        {
#ifdef VOLESTI_DEBUG
            std::cout<<"Could not construct constaints for the Linear Program for membership "<<e<<std::endl;
#endif
            return false;
        }
    }

    set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */

    for(j=0; j<k; j++){
        //colno[j] = j+1; // last column
        mem_comb[j] = 0.0;
        set_bounds(lp, j+1, -infinite, infinite);
    }

    /* add the row to lpsolve */
    try {
        if(!set_obj_fnex(lp, k, mem_comb, colno)) throw false;
    }
    catch (bool e)
    {
#ifdef VOLESTI_DEBUG
        std::cout<<"Could not construct objective function for the Linear Program for membership "<<e<<std::endl;
#endif
        return false;
    }

    //set the bounds

    /* set the object direction to maximize */
    set_maxim(lp);

    /* I only want to see important messages on screen while solving */
    set_verbose(lp, NEUTRAL);

    /* Now let lpsolve calculate a solution */
    if (solve(lp) != OPTIMAL){
        delete_lp(lp);
        return false;
    }
    get_variables(lp, mem_comb);
    for (int l = 0; l < k; ++l) {
        //std::cout<<"row["<<l<<"] = "<<mem_comb[l]<<std::endl;
    }
    //std::cout<<"\n";
    delete_lp(lp);
    return true;

}



template <class MT, class VT, class Point, typename NT2>
void get_feas_point(MT &T, MT &A, VT &b, MT &Aeq, VT &beq, Point &q, NT2 *mem_comb) {

    typedef typename Point::FT NT;
    int d=q.dimension(), m = Aeq.rows(), k = T.cols(), m2 = A.rows();
    lprec *lp;
    int Ncol=k, *colno = NULL, j, i, Nrows = d+m+m2;

    //REAL *row = NULL;

    try
    {
        lp = make_lp(Nrows, Ncol);
        if(lp == NULL) throw false;
    }
    catch (bool e) {
#ifdef VOLESTI_DEBUG
        std::cout<<"Could not construct Linear Program for membership "<<e<<std::endl;
#endif
        return;
    }

    REAL infinite = get_infinite(lp); /* will return 1.0e30 */

    try
    {
        colno = (int *) malloc(Ncol * sizeof(*colno));
        //row = (REAL *) malloc(Ncol * sizeof(*row));
    }
    catch (std::exception &e)
    {
#ifdef VOLESTI_DEBUG
        std::cout<<"Linear Program for membership failed "<<e.what()<<std::endl;
#endif
        return;
    }

    set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */

    for (i = 0;  i< d; ++i) {
        /* construct all rows */
        for(j=0; j<k; j++){
            colno[j] = j+1;
            mem_comb[j] = T(i,j);
        }

        /* add the row to lpsolve */
        try {
            if(!add_constraintex(lp, k, mem_comb, colno, EQ, q[i])) throw false;
        }
        catch (bool e)
        {
#ifdef VOLESTI_DEBUG
            std::cout<<"Could not construct constaints for the Linear Program for membership "<<e<<std::endl;
#endif
            return;
        }
    }

    for (i = 0;  i< m; ++i) {
        /* construct all rows */
        for(j=0; j<k; j++){
            //colno[j] = j+1;
            mem_comb[j] = Aeq(i,j);
        }

        /* add the row to lpsolve */
        try {
            if(!add_constraintex(lp, k, mem_comb, colno, EQ, beq(i))) throw false;
        }
        catch (bool e)
        {
#ifdef VOLESTI_DEBUG
            std::cout<<"Could not construct constaints for the Linear Program for membership "<<e<<std::endl;
#endif
            return;
        }
    }

    for (i = 0;  i< m2; ++i) {
        /* construct all rows */
        for(j=0; j<k; j++){
            //colno[j] = j+1;
            mem_comb[j] = A(i,j);
        }

        /* add the row to lpsolve */
        try {
            if(!add_constraintex(lp, k, mem_comb, colno, LE, b(i))) throw false;
        }
        catch (bool e)
        {
#ifdef VOLESTI_DEBUG
            std::cout<<"Could not construct constaints for the Linear Program for membership "<<e<<std::endl;
#endif
            return;
        }
    }

    set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */

    for(j=0; j<k; j++){
        //colno[j] = j+1; // last column
        mem_comb[j] = 0.0;
        set_bounds(lp, j+1, -infinite, infinite);
    }

    /* add the row to lpsolve */
    try {
        if(!set_obj_fnex(lp, k, mem_comb, colno)) throw false;
    }
    catch (bool e)
    {
#ifdef VOLESTI_DEBUG
        std::cout<<"Could not construct objective function for the Linear Program for membership "<<e<<std::endl;
#endif
        return;
    }

    //set the bounds

    /* set the object direction to maximize */
    set_maxim(lp);

    /* I only want to see important messages on screen while solving */
    set_verbose(lp, NEUTRAL);

    /* Now let lpsolve calculate a solution */
    if (solve(lp) != OPTIMAL){
        std::cout<<"not feasible"<<std::endl;
        delete_lp(lp);
        return;
    }
    get_variables(lp, mem_comb);
    for (int l = 0; l < k; ++l) {
        std::cout<<"feas_point_row["<<l<<"] = "<<mem_comb[l]<<std::endl;
    }
    std::cout<<"\n";
    delete_lp(lp);
    return;

}



#endif
