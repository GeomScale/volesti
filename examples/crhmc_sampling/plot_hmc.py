# VolEsti (volume computation and sampling library)
#
# Copyright (c) 2012-2020 Vissarion Fisikopoulos
# Copyright (c) 2018-2020 Apostolos Chalkis
# Copyright (c) 2020-2020 Marios Papachristou
#
# Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.
#
# Licensed under GNU LGPL.3, see LICENCE file

import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns
import math
sns.set_style("white")
from io import StringIO

if __name__ == '__main__':

    usage = '''
        Script to plot samples from logconcave density
        Example usage:

        ./simple_hmc | python3 plot_hmc.py

    '''

    argparser = argparse.ArgumentParser(usage=usage)
    argparser.add_argument('--name', type=str, default='hmc_plot', help='Plot name')
    argparser.add_argument('--limits', type=str, default='auto', help='Plot limits for scatter plots')
    argparser.add_argument('--save', action='store_true', help='Save output figures')
    argparser.add_argument('--max_marginals', default=-1, type=int, help='Plot the maximum number of marginals')
    args = argparser.parse_args()

    temp = StringIO(sys.stdin.read())
    data = np.loadtxt(temp)

    if len(data.shape) != 2:
        data = data.reshape(data.shape[-1], 1)

    if args.max_marginals > 0:
        dims = min(data.shape[-1], args.max_marginals)
    else:
        dims = data.shape[-1]

    print('Number of dimensions: {}'.format(dims))

    print('Plotting histograms of marginal densities')
    fig, ax2d = plt.subplots(ncols=dims, nrows=1, squeeze=False)
    axli = ax2d.flatten()

    plt.suptitle('{}: Marginals'.format(args.name))
    #t = np.linspace(-1, 1, 1000)
    #a=(1/1.31649)*np.exp(-2*t*t-t)
    for i in range(1, 1 + dims):
        sns.histplot(data[:, i-1], bins=50, color="orange", ax=axli[i-1], stat='probability', kde=True)
        #axli[i-1].plot(t,a,'r')
        axli[i-1].set_xlabel('$x_{}$'.format(i))
        axli[i-1].set_ylabel('$\pi(x_{})$'.format(i))

    if args.save:
        plt.savefig('{}_marginals.png'.format(args.name))

    print('Plotting samples')
    fig, ax2d = plt.subplots(ncols=dims, nrows=1, squeeze=False)
    axli = ax2d.flatten()

    plt.suptitle('{}: Samples'.format(args.name))
    for i in range(1, 1 + dims):
        axli[i-1].plot(data[:, i-1])
        axli[i-1].set_xlabel('Number of samples (t)')
        axli[i-1].set_ylabel('$x_{}(t)$'.format(i))

    if args.save:
        plt.savefig('{}_samples.png'.format(args.name))

    if dims in [2, 3]:
        print('Plotting scatterplot')
        fig = plt.figure()
        if dims == 2:
            ax = fig.add_subplot(111)
            ax.scatter(data[:, 0], data[:, 1])
        if dims == 3:
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(data[:, 0], data[:, 1], data[:, 2])

        ax.set_xlabel('$x_1$')
        ax.set_ylabel('$x_2$')
        if dims == 3:
            ax.set_zlabel('$x_3$')

        fig.suptitle('{}: {}D Scatter Plot'.format(args.name, dims))

        if args.limits != 'auto':
            try:
                limits = int(args.limits)
                ax.set_xlim(-limits, limits)
                ax.set_ylim(-limits, limits)
                if dims == 3:
                    ax.set_zlim(-limits, limits)
            except ValueError:
                pass

        if args.save:
            plt.savefig('{}_scatter.png'.format(args.name))

    plt.show()
