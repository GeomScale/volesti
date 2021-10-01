'''
    VolEsti (volume computation and sampling library)

    Copyright (c) 2012-2021 Vissarion Fisikopoulos
    Copyright (c) 2020-2021 Apostolos Chalkis
    Copyright (c) 2021- Marios Papachristou

    Script to convert post processed .mat file to .ine format
    in order to be processed by the C++ backend.

    The polytope is given in a .mat file and is represented in
    the form Ax <= b. Moreover, the Chebyshev ball details (center, radius)
    are given.

'''
import argparse
import numpy as np
import scipy.io
import os
import sys

def get_args():
    parser = argparse.ArgumentParser(description='Converts a polytope given in canonical form (post-processed) in .mat format  to an .ine file')
    parser.add_argument('-i', type=str, help='Input .mat file. Must contain entries A, b, center, and radius.')
    parser.add_arugment('-o', default=None, help='Output file prefix. Same as input file (with .ine extension) if left unspecified')

    return parser.parse_args()

if __name__ == '__main__':
    args = get_args()
    filename = args.i
    outfilename = args.o
    mat = scipy.io.loadmat(filename)

    try:
        A = mat['polytope']['A'][0, 0]
        b = mat['polytope']['b'][0, 0]
        center = mat['polytope']['center'][0, 0]
        radius = mat['polytope']['radius'][0, 0]
    except KeyError:
        raise Exception('Uncompatible .mat file')
        sys.exit(1)

    b = b - A @ center

    C = np.hstack((b, - A))

    m, d = C.shape

    if args.o:
        outfile = args.o + '.ine'
        outfile_inner_ball = args.o + '.inner_ball'
    else:
        outfile = os.path.splitext(filename)[0] + '.ine'
        outfile_inner_ball = os.path.splitext(filename)[0] + '.inner_ball'

    with open(outfile, 'w+') as f:
        f.write('filename\n')
        f.write('begin\n')
        f.write('{} {} real\n'.format(m, d))
        np.savetxt(f, C, fmt="%.8f")
        f.write('\n')
        f.write('end\n')
        f.write('input_incidence\n')

    with open(outfile_inner_ball, 'w+') as f:
        np.savetxt(f, center, fmt="%.10f")
        np.savetxt(f, radius, fmt="%.10f")
