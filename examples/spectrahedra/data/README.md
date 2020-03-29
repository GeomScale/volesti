## SDPA format

The file sdp_n2m3.txt contains a semidefinite program in SDPA format. For details in the SDPA format you can
search the official manual.

The contents of the file are:

2 <br>
1<br>
3<br>
1 1<br>
-1  0  0<br>
 0 -2  1<br>
 0  1 -2<br>
 1 -0 -0<br>
-0 -0 -1<br>
-0 -1 -0<br>
-0 -0  1<br>
-0 -0 -0<br>
 1 -0 -0
 
It represents a spectrahedron in 2 dimensions and a linear objective function. The spectrahedron is described
by a linear matrix inequality, of the form x1 A1 + ... + xn An - A0 >= 0, where >= 0 denotes positive semidefiniteness.

### Explanation line by line
- 2 : The number of dimensions. The vector of the objective function is of length 2 and we need 3 matrices for 
   the linear matrix inequality, i.e. x A1 + y A2 - A0 >= 0.
- 1 : Please lookup SDPA format for more details. Currently support is only for value 1.
- 3 : The size of the matrices, i.e. 3x3
- 1 1 : The vector of the objective function, i.e. (1,1) or f(x,y) = x + y
- -1  0  0: The first row of A0
- 0 -2  1: The second row of A0
- 0  1 -2: The third row of A0
- 1 -0 -0: The first row of A1
- and so on, till all 3 matrices are defined