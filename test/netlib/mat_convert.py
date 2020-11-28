import scipy.io
import os
import sys
import numpy as np

filename = sys.argv[-1]

mat = scipy.io.loadmat(filename)
eps = 1e-2
val = 9999999

Aeq = mat['problem']['Aeq'][0, 0].toarray()
beq = mat['problem']['beq'][0, 0]
lb = mat['problem']['lb'][0, 0]
ub = mat['problem']['ub'][0, 0]

lb[lb == - np.inf] = -val
lb[lb == np.inf] = val
ub[ub == - np.inf] = -val
ub[ub == np.inf] = val

# Equality constraints Ax = b
# are replaced by Ax <= b and - Ax <= - b
# Lower bounds are defined by -Ix <= -lb
# .ine files expect a format (b_hat, -A_hat)
C = np.vstack (
    (
        np.hstack((eps + beq, -Aeq)),
        np.hstack((eps -beq, Aeq)),
        np.hstack((-lb, np.eye(lb.shape[0]))),
        np.hstack((ub, -np.eye(ub.shape[0])))
    )
)

m, d = C.shape

outfile = os.path.splitext(filename)[0] + '.ine'

with open(outfile, 'w+') as f:
    f.write('filename\n')
    f.write('begin\n')
    f.write('{} {} real\n'.format(m, d))
    np.savetxt(f, C, fmt="%.4f")
    f.write('\n')
    f.write('end\n')
    f.write('input_incidence\n')
