import math
import numpy as np
from matplotlib import pyplot as plt
import random

d=2
n=2000
mu=np.array([-0.5,0.5]).reshape(d,1)
# the mean of the two gaussian is -0.5 and 0.5
y=np.concatenate([np.random.normal(mu[i],0.1,n) for i in range(d)])
# standard deviation is 0.1
plt.hist(y,bins=100)
for i in y:
    print(str(i))
plt.show()