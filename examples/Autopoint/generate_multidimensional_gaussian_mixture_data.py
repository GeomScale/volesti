import numpy as np
import random
d=4 # dimension of the data
n=4000
mean1_=[random.random() for i in range(d) ]
#the two mean are directly opposite for maximum sampling efficiency
mean2_=[-i for i in mean1_]
cov = np.identity(d)
pts1 = np.random.multivariate_normal(mean1_, cov, size=n)
pts2 = np.random.multivariate_normal(mean2_, cov, size=n)
sample=np.r_[pts1,pts2]
for i in sample:
    print(" ".join(map(str, i)))
import matplotlib.pyplot as plt
plt.plot(sample[:, 0], sample[:, 1], '.', alpha=0.5)
plt.axis('equal')
plt.grid()
plt.show()

