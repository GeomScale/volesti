Practical volume computation of structured convex bodies, and an application to modeling
                      portfolio dependencies and financial crises.
   Ludovic Cales, Apostolos Chalkis, Ioannis Z. Emiris and Vissarion Fisikopoulos

By Tolis Chalkis, N.K. Univ. Athens, Greece, 2017.

Here we give all the implementations that are mentioned in the paper.

------------------ 
1. Compile sources 
------------------ 
 
In folder project_main execute: 
 
$ cgal_create_CMakeLists -s executable 
$ cmake -DCGAL_DIR=_YOUR_CGAL_PATH_ . 
$ make 

where _YOUR_CGAL_PATH_ is the path where CGAL library was compiled.

We are using CGAL 11 and gcc6.

-------------------------
2. How to run algorithms
-------------------------

./exe -h   would give you some details about the algorithms you are able to run and how to call them.

Examples:

i) ./exe -r 1 [dim] [k]  
will compare rejection for arbitrary and the unit simplex for dimension equals to dim and number of points sampling equals to k


ii) ./exe -r 2 [dim] [k] [file.txt]   will run rejection for two hypeplanes intersecting the unit simplex. The file.txt must contain the two hyperplanes in two rows. In every row coofficients go first and the constant zi last. For example row:

2 -4 7 12 

is the 3-dimensional hyperplane: 2x1 -4x2 +7x3 <= 12

"dim29_two_hyperplanes_for_unit.txt" is such a txt file.


iii) ./exe -r 3 [dim] [k] [file.txt]     
Will run rejection for two families of parallel hyperplanes intersecting the facet simplex. The file.txt must contain two rows. Each row must contain the coefficients that form the direction of each family. "dim100_two_families_par_hyp.txt" is such a txt file.


iv) ./exe -r 4 [dim] [k] [file1.txt] [file2.txt] 
Will run rejection for a family of parallel hyperplanes and a family of cocentric ellipsoids intersecting the facet simplex and centred at the origin. The file1.txt must contain the family of parallel hyperplanes as "dim100_one_fam_hyp.txt". The file2.txt must contain the matrix of the ellipsoids. "dim100_mat_ellips.txt" is such a txt file.


v) ./exe -ve 1 1 [dim] [k] [W] [e] [file.txt]
Will run VolEsti for an ellipsoid intersecting the unit simplex starting with an almost random inscribed ball. The file.txt must contain the matrix of the ellipsoid and the last dim+1 numbers correspond to the center and the constant of the ellipsoid. "dim10_mat_ellips_ve.txt" is such a file. Rejection will run as well.


vi) ./exe -ve 1 2 [dim] [k] [W] [e] [file1.txt] [file2.txt] 
Will run VolEsti for an ellipsoid intersecting the unit simplex starting with inscribed ball from conical optimization problem. The file1.txt must contain the matrix of the ellipsoid and the last dim+1 numbers correspond to the center and the constant of the ellipsoid. "dim10_mat_ellips_ve.txt" is such a file. The file2.txt must contain dim+1 numbers corresponding to the center and the radius of the inscribed ball, "dim10_center_radius_ve.txt" is such a txt. Rejection will run as well.
To compute the center and the radius of the the inscribed ball you have to use matlab file in matlab folder. You have to give as inputs the dimension and the number of observations and an ellipsoid will be constructed from real data centred at the center of the unit simplex.


vii) ./exe -ve 2 [dim] [k] [W] [e] [file1.txt] [file2.txt] 
Will run VolEsti for an ellipsoid and two parallel hyperplanes intersecting the unit simplex. The file1.txt must contain the ellipsoid matrix, the center and the constant. The file2.txt must contain one row. The dim first numbers are the coefficients of the hyperplanes and the last two are the constants. "dim10_two_hyp_par.txt" is such a txt. Rejection will run as well.

viii) ./exe -ve 3 [dim] [k] [W] [e] [file1.txt] [file2.txt]
Will run VolEsti for two parallel hyperplanes and two cocentric ellipsoids intersescting the unit simplex which is a non convex body. The file1.txt must contain the matrix of the ellipsoids, the center and the two last numbers correspond to the two constants c1>c2. "dim10_two_ellips_nonconv.txt" is such a txt. The file2.txt must contain the coeefixcients and the constants of the two parallel hyperplanes. Rejection will run as well.


ix) ./exe -varsi [dim] [file.txt] 
Will run Varsi formula. The file.txt must contain a hyperplane which intersects with the unit simplex. The "dim10_hyperplane.txt" is such a txt.


x) ./exe -Lawn 1 [dim] [file.txt]
Will run Lawrence formula, with multiprecision computations, for two parallel hyperplanes intersecting the unbit simplex. Rejection will run as well.


xi) ./exe -Lawn 2 [dim] [file.txt]
Same with (x) but with double precision computations.
