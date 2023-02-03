import numpy as np
import random
import argparse

argParser = argparse.ArgumentParser()
argParser.add_argument("-d", "--dimension", type=int,help="dimension of mean of normal distribution.There will be samples from two normal distribution")
argParser.add_argument("-n", "--number", type=int,default=4000,help="number of samples")
args = argParser.parse_args()
mean1_=[random.random() for i in range(args.dimension) ]
#the two mean are directly opposite for maximum sampling efficiency
mean2_=[-i for i in mean1_]
cov = 0.1*np.identity(args.dimension)
pts1 = np.random.multivariate_normal(mean1_, cov, size=args.number//2)
pts2 = np.random.multivariate_normal(mean2_, cov, size=args.number//2)
sample=np.r_[pts1,pts2]
for i in sample:
    print(" ".join(map(str, i)))

