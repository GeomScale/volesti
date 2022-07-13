#!/usr/bin/bash
LPSOLVE='/usr/lib/lp_solve/liblpsolve55.so'

cmake . -DLP_SOLVE=$LPSOLVE 

make

./sampler

# cat uniform_sampling.txt
