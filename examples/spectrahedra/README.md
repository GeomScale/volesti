# Examples for Spectrahedra

## Compilation
In folder examples, first run cmake, to create the makefile:

```bash
cmake .
```

Then, in folder examples/spectrahedra compile using the makefile:

```bash
make
```

## List of examples
- Example 1: Read/write SDPA format files

## Examples
### Example 1: Read/write SDPA format files

In this example, we will read a semidefinite program from a SDPA format input file, print it
and then write it to a new SDPA format file. Run the example with:

```bash
./readWriteSdpaFile
```

The input file is data/sdp_n2m3.txt. It contains a semidefinite program in SDPA format. A semidefinite program
(dual form) consists of a linear matrix inequality (describing a spectrahedron) of the form A1 + ... + xn An - A0 >= 0,
where >= 0 denotes positive semidefiniteness. For details in the SDPA format you can search the official manual.

The contents of the file are:

```bash
2
1
3
1 1
-1  0  0
 0 -2  1
 0  1 -2
 1 -0 -0
-0 -0 -1
-0 -1 -0
-0 -0  1
-0 -0 -0
 1 -0 -0
```

 
It represents a spectrahedron in 2 dimensions, described by a linear matrix inequality with
3x3 matrices, and a linear objective function. 

##### Explanation line by line
- 2 : The number of dimensions. The vector of the objective function is of length 2 and we need 3 matrices for 
   the linear matrix inequality, i.e. x A1 + y A2 - A0 >= 0.
- 1 : Please lookup SDPA format for more details. Currently, support is only for value 1.
- 3 : The size of the matrices, i.e. 3x3
- 1 1 : The vector of the objective function, i.e. (1,1) or f(x,y) = x + y
- -1  0  0: The first row of A0
- 0 -2  1: The second row of A0
- 0  1 -2: The third row of A0
- 1 -0 -0: The first row of A1
- and so on, till all 3 matrices are defined