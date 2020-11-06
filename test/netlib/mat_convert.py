import scipy.io
import os
import sys
import numpy as np

filename = sys.argv[-1]

mat = scipy.io.loadmat(filename)

A = - mat['problem']['Aeq'][0, 0].toarray()
b = mat['problem']['beq'][0, 0]

C = np.hstack((b, A))

m, d = A.shape

print('filename')
print('begin')
print('{} {} real'.format(m, d + 1))
np.savetxt(sys.stdout, C, fmt="%.4f")

print('end')
print('input_incidence')

