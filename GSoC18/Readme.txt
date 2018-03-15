Sampling and volume approximation, project for GSoC 2018.

By Tolis Chalkis, N.K. Univ. Athens, Greece, 2018.

Here we give the three answers for the tests are required.

------------------ 
1. Easy: use VolEsti to compute the volume of the 10-dimensional hypercube.
------------------ 

We compiled the project in folder examples and use the file "cube10.ine" which contains the the 10-dimensional hypercube in H-represantation.
We know that the exact volume of the corresponding hypercube is 2^10=1024.

We have excecuted the command below:

./vol -f1 cube10.ine
  Reading input from file...
  Experiment 1
  10 20 1 -1 1 [-0,-2] 9210 11 990.7898 [990.7898,990.7898] 0 991.7898 0 0.279661 0.001171 

The computed volume is 990.7898.

VolEsti uses the default error and walk length. So we could give different values as below:

./vol -f1 cube10.ine -e 0.2 -w 12
  Reading input from file...
  Experiment 1 
  10 20 1 -1 0.2 [-0.8,-1.2] 230258 12 1020.839 [1020.839,1020.839] 0 1021.839 0 7.144641 0.000354 

The computed volume is 1020.839 which is closer to the exact volume as we set a larger number of points in the Hit-And-Run in each step and a larger number for the walk length.

Moreover we have the option to declare that the input polytope is a hypercube:

./vol --cube -f1 cube10.ine -e 0.3 -w 11
  Reading input from file...
  Experiment 1 
  10 20 1 1 0.3 [0.7,1.3] 102337 11 1015.48 [1015.48,1015.48] 0 -1014.48 0 2.962104 0.001155

The computed volume is 1015.48

------------------
2. Medium: write a function that samples using the ball walk (see page 3 of https://www.cc.gatech.edu/~vempala/papers/survey.pdf). 
Modify VolEsti to compute volumes using your function.
------------------

In vol.cpp we add a new option for ball walk. If you run the command:

./vol -f1 cube10.ine -bw

Then ball walk is going to be used.

In random_samplers.h we implement a function for ball walk. Moreover in vol_rand.h we add to the structure var two more members: delta and ball_walk. The first is a double for the radius that ball walk algorithm uses and the second is a boolean variable that is true if -bw is given as an input string, otherwise is false.
