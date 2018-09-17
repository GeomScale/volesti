#install.packages("Peacock.test")

library(volesti)
library(Peacock.test)

#--------------------------------------------------------------------#
#-------------Uniformity test in 3d unit simplex---------------------#
#--------------------------------------------------------------------#
print('- - - - - - - - - - -  - - - - - - - - - - -- - - - - - - - -')
print('Uniformity test in 3d unit simplex')
print('- - - - - - - - - - - - - -- - - - - - - - - - -  - - - - - -')
print(' ')


# Generate the 3d unit simplex
PolyList = GenSimplex(3, 'H')
# sample 1000 uniform points from the 3d unit simplex
sample1 = sample_simplex(dimension = 3, N = 1000)


#-----Coordinate Directions Hit-and-run-------#

# sample 1000 points using coordinate directions hit-and-run with walk step equals to 1
sample2 = sample_points(A = PolyList$A, b = PolyList$b, N = 1000, walk_length = 1)
print(paste0('[CDHR][W=',1,'] the value of the test statistic = ',peacock3(sample1, sample2)))

# sample 1000 points using coordinate directions hit-and-run with walk step equals to 10
sample2 = sample_points(A = PolyList$A, b = PolyList$b, N = 1000, walk_length = 10)
print(paste0('[CDHR][W=',10,'] the value of the test statistic = ',peacock3(sample1, sample2)))
print('- - - - - - - - - - - - - - - - - - - -')
#---------------------------------------------------------------------------#

#-----Random Directions Hit-and-run-------#

# sample 1000 points using random directions hit-and-run with walk step equals to 1
sample2 = sample_points(A = PolyList$A, b = PolyList$b, N = 1000, walk_length = 1, coordinate = FALSE)
print(paste0('[RDHR][W=',1,'] the value of the test statistic = ',peacock3(sample1, sample2)))

# sample 1000 points using random directions hit-and-run with walk step equals to 10
sample2 = sample_points(A = PolyList$A, b = PolyList$b, N = 1000, walk_length = 10, coordinate = FALSE)
print(paste0('[RDHR][W=',10,'] the value of the test statistic = ',peacock3(sample1, sample2)))
print('- - - - - - - - - - - - - - - - - - - -')
#---------------------------------------------------------------------------#

#-----Ball walk-------#

# sample 1000 points using ball walk with walk step equals to 1
sample2 = sample_points(A = PolyList$A, b = PolyList$b, N = 1000, walk_length = 1, ball_walk = TRUE)
print(paste0('[BW][W=',1,'] the value of the test statistic = ',peacock3(sample1, sample2)))

# sample 1000 points using ball walk with walk step equals to 10
sample3 = sample_points(A = PolyList$A, b = PolyList$b, N = 1000, walk_length = 10, ball_walk = TRUE)
print(paste0('[BW][W=',10,'] the value of the test statistic = ',peacock3(sample1, sample2)))
print('- - - - - - - - - - - - - - - - - - - -')
print(' ')
#---------------------------------------------------------------------------#


#--------------------------------------------------------------------#
#-------------Uniformity test in 2d skinny simplex-------------------#
#--------------------------------------------------------------------#
print('- - - - - - - - - - -  - - - - - - - - - - -- - - - - - - - -')
print('Uniformity test in a 2d skinny simplex')
print('- - - - - - - - - - - - - -- - - - - - - - - - -  - - - - - -')
print(' ')


# V-representation of simplex
V  = matrix(c(0,0,0,7,100,0), ncol=2, nrow=3, byrow=TRUE)
# H-representation of simplex
A = matrix(c(-1,0,0,-1,0.07,1), ncol=2, nrow=3, byrow=TRUE)
b=c(0,0,7)

sample1 = sample_simplex(vertices =V, N = 100)
sample2 = sample_points(A=A, b=b, N=100, walk_length = 1)
print(paste0('[CDHR][W=',1,'] the value of the test statistic = ',peacock2(sample1, sample2)))
print(' ')

sample2 = sample_points(A=A, b=b, N=100, walk_length = 10)
print(paste0('[CDHR][W=',10,'] the value of the test statistic = ',peacock2(sample1, sample2)))
print(' ')

sample2 = sample_points(A=A, b=b, N=100, walk_length = 50)
print(paste0('[CDHR][W=',50,'] the value of the test statistic = ',peacock2(sample1, sample2)))
print(' ')
print(' - - - - - - - - - ')

sample2 = sample_points(A=A, b=b, N=100, walk_length = 1, coordinate = FALSE)
print(paste0('[RDHR][W=',1,'] the value of the test statistic = ',peacock2(sample1, sample2)))
print(' ')

sample2 = sample_points(A=A, b=b, N=100, walk_length = 10, coordinate = FALSE)
print(paste0('[RDHR][W=',10,'] the value of the test statistic = ',peacock2(sample1, sample2)))
print(' ')

sample2 = sample_points(A=A, b=b, N=100, walk_length = 50, coordinate = FALSE)
print(paste0('[RDHR][W=',50,'] the value of the test statistic = ',peacock2(sample1, sample2)))
print(' ')
print(' - - - - - - - - - ')

sample2 = sample_points(A=A, b=b, N=100, walk_length = 1, ball_walk = TRUE)
print(paste0('[BW][W=',1,'] the value of the test statistic = ',peacock2(sample1, sample2)))
print(' ')

sample2 = sample_points(A=A, b=b, N=100, walk_length = 10, ball_walk = TRUE)
print(paste0('[BW][W=',10,'] the value of the test statistic = ',peacock2(sample1, sample2)))
print(' ')

sample2 = sample_points(A=A, b=b, N=100, walk_length = 50, ball_walk = TRUE)
print(paste0('[BW][W=',50,'] the value of the test statistic = ',peacock2(sample1, sample2)))
print(' ')
print(' - - - - - - - - - ')

#--------------------------------------------------------------------------#
#---------------Coordinate Dir HnR -vs- Random Dir HnR---------------------#
#--------------------------------------------------------------------------#
print('- - - - - - - - - - -  - - - - - - - - - - -- - - - - - - - -')
print('Compare Coordinate Directions HnR with Random Directions HnR')
print('- - - - - - - - - - - - - -- - - - - - - - - - -  - - - - - -')
print(' ')

# Generate a 3d cube
print('--3d cube--')
PolyList = GenCube(3, 'H')
# sample 5e03 points using Coordinate Dir HnR
sample1 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, verbose = TRUE)
# sample 5e03 points using Random Dir HnR
sample2 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, coordinate = FALSE, verbose = TRUE)
print('- - - - - - - - - - - - - - - - - - - -')

# Generate a 3d cross polytope
print('--3d cross polytope--')
PolyList = GenCross(3, 'H')
# sample 5e03 points using Coordinate Dir HnR
sample1 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, verbose = TRUE)
# sample 5e03 points using Random Dir HnR
sample2 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, coordinate = FALSE, verbose = TRUE)
print('- - - - - - - - - - - - - - - - - - - -')

# Generate a 10d cube
print('--10d cube--')
PolyList = GenCube(10, 'H')
# sample 5e03 points using Coordinate Dir HnR
sample1 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, verbose = TRUE)
# sample 5e03 points using Random Dir HnR
sample2 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, coordinate = FALSE, verbose = TRUE)
print('- - - - - - - - - - - - - - - - - - - -')

# Generate a 10d cross polytope
print('--10d cross polytope--')
PolyList = GenCross(10, 'H')
# sample 5e03 points using Coordinate Dir HnR
sample1 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, verbose = TRUE)
# sample 5e03 points using Random Dir HnR
sample2 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, coordinate = FALSE, verbose = TRUE)
print('- - - - - - - - - - - - - - - - - - - -')

# Generate a 30d cube
print('--30d cube--')
PolyList = GenCube(30, 'H')
# sample 5e03 points using Coordinate Dir HnR
sample1 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, verbose = TRUE)
# sample 5e03 points using Random Dir HnR
sample2 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, coordinate = FALSE, verbose = TRUE)
print('- - - - - - - - - - - - - - - - - - - -')

# Generate a 15d cross polytope
print('--15d cross polytope--')
PolyList = GenCross(15, 'H')
# sample 5e03 points using Coordinate Dir HnR
sample1 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, verbose = TRUE)
# sample 5e03 points using Random Dir HnR
sample2 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, coordinate = FALSE, verbose = TRUE)
print('- - - - - - - - - - - - - - - - - - - -')

# Generate a 80d cube
print('--80d cube--')
PolyList = GenCube(80, 'H')
# sample 5e03 points using Coordinate Dir HnR
sample1 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, verbose = TRUE)
# sample 5e03 points using Random Dir HnR
sample2 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, coordinate = FALSE, verbose = TRUE)
print('- - - - - - - - - - - - - - - - - - - -')

# Generate a 10d cross polytope
print('--10d cross polytope--')
PolyList = GenCross(10, 'H')
# sample 5e03 points using Coordinate Dir HnR
sample1 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, verbose = TRUE)
# sample 5e03 points using Random Dir HnR
sample2 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, coordinate = FALSE, verbose = TRUE)
print('- - - - - - - - - - - - - - - - - - - -')


#--------------------------------------------------------------------------#
#---------------Coordinate Dir HnR -vs- Ball walk---------------------#
#--------------------------------------------------------------------------#
print('- - - - - - - - - - -  - - - - - - - - - - -- - - - - - - - -')
print('Compare Coordinate Directions HnR with Ball walk')
print('- - - - - - - - - - - - - -- - - - - - - - - - -  - - - - - -')
print(' ')

# Generate a 3d cube
print('--3d cube--')
PolyList = GenCube(3, 'H')
# sample 5e03 points using Coordinate Dir HnR
sample1 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, verbose = TRUE)
# sample 5e03 points using Random Dir HnR
sample2 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, ball_walk = TRUE, verbose = TRUE)
print('- - - - - - - - - - - - - - - - - - - -')

# Generate a 3d cross polytope
print('--3d cross polytope--')
PolyList = GenCross(3, 'H')
# sample 5e03 points using Coordinate Dir HnR
sample1 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, verbose = TRUE)
# sample 5e03 points using Random Dir HnR
sample2 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, ball_walk = TRUE, verbose = TRUE)
print('- - - - - - - - - - - - - - - - - - - -')

# Generate a 10d cube
print('--10d cube--')
PolyList = GenCube(10, 'H')
# sample 5e03 points using Coordinate Dir HnR
sample1 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, verbose = TRUE)
# sample 5e03 points using Random Dir HnR
sample2 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, ball_walk = TRUE, verbose = TRUE)
print('- - - - - - - - - - - - - - - - - - - -')

# Generate a 10d cross polytope
print('--10d cross polytope--')
PolyList = GenCross(10, 'H')
# sample 5e03 points using Coordinate Dir HnR
sample1 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, verbose = TRUE)
# sample 5e03 points using Random Dir HnR
sample2 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, ball_walk = TRUE, verbose = TRUE)
print('- - - - - - - - - - - - - - - - - - - -')

# Generate a 30d cube
print('--30d cube--')
PolyList = GenCube(30, 'H')
# sample 5e03 points using Coordinate Dir HnR
sample1 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, verbose = TRUE)
# sample 5e03 points using Random Dir HnR
sample2 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, ball_walk = TRUE, verbose = TRUE)
print('- - - - - - - - - - - - - - - - - - - -')

# Generate a 15d cross polytope
print('--15d cross polytope--')
PolyList = GenCross(15, 'H')
# sample 5e03 points using Coordinate Dir HnR
sample1 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, verbose = TRUE)
# sample 5e03 points using Random Dir HnR
sample2 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, ball_walk = TRUE, verbose = TRUE)
print('- - - - - - - - - - - - - - - - - - - -')

# Generate a 80d cube
print('--80d cube--')
PolyList = GenCube(80, 'H')
# sample 5e03 points using Coordinate Dir HnR
sample1 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, verbose = TRUE)
# sample 5e03 points using Random Dir HnR
sample2 = sample_points(A = PolyList$A, b = PolyList$b, N = 5000, ball_walk = TRUE, verbose = TRUE)
print('- - - - - - - - - - - - - - - - - - - -')



