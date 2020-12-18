import scipy.io
import os
import sys
import numpy as np

filename = sys.argv[-1]

form = 'canonical'

mat = scipy.io.loadmat(filename)
eps = 1e-2
scale = 2
val = sys.maxsize

# Distinguish between inputs of type (b, -A) (sign = -1) and (b, A) (sign = 1)
sign = -1

if form == 'canonical':
    # Polytope is provided in the form Ax <= b and x >= 0
    A = mat['polytope']['A'][0, 0]
    b = mat['polytope']['b'][0, 0]
    center = mat['polytope']['center'][0, 0]
    radius = mat['polytope']['radius'][0, 0]

    C = np.hstack((b + 1, sign * A))

if form == 'standard':
    # Polytope is provided in the form A_eq x = b_eq and x >= lb and x <= ub
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
            np.hstack((eps + beq, sign * Aeq)),
            np.hstack((eps -beq, - sign * Aeq)),
            np.hstack((-lb, - sign * np.eye(lb.shape[0]))),
            np.hstack((ub, sign * np.eye(ub.shape[0])))
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

outfile_inner_ball = os.path.splitext(filename)[0] + '.inner_ball'

with open(outfile_inner_ball, 'w+') as f:
    np.savetxt(f, center, fmt="%.4f")
    np.savetxt(f, radius, fmt="%.4f")
