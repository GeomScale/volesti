---
title: 'volesti: C++ library for sampling and volume computation on convex bodies'
tags:
  - C++
  - geometry
  - randomization
  - Monte-Carlo methods
  - convexity
authors:
  - name: Apostolos Chalkis
    orcid: 0000-0000-0000-0000
    equal-contrib: true
    affiliation: "2, 4" # (Multiple affiliations must be quoted)
  - name: Vissarion Fisikopoulos
    corresponding: true # (This is how to denote the corresponding author)
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: "1, 4"
  - name: Marios Papachristou
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 5
  - name: Elias Tsigaridas
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: "3, 4"
affiliations:
 - name: National & Kapodistrian University of Athens, Greece
   index: 1
 - name: Quantagonia
   index: 2
 - name: Inria Paris and IMJ-PRG, Sorbonne Universit\`e and Paris Universit\`e
   index: 3
 - name: GeomScale
   index: 4
 - name: Cornell University
   index: 5
date: 11 March 2024
bibliography: paper.bib

---

# Summary

Sampling from high-dimensional distributions and volume approximation of convex
bodies are fundamental operations that appear in optimization, finance,
engineering, artificial intelligence, and machine learning.
In this paper, we present `volesti`, a C++ library that provides efficient, randomized algorithms for
volume estimation---a special case of integration---as well as general logconcave sampling and
rounding for convex bodies.
Since the implemented methods are Monte-Carlo algorithms the library also provides several
high-dimensional statistical tests.

`volesti`'s focus is scalability in high dimensions that could be in the order of hundreds or thousands
depending on the problem.
Another novelty of the library is the variety of handling inputs.
`volesti` supports three different types of polyhedra [@Ziegler:1995], specrahedra [Ramana:1999]
and general non-linear convex objects.

`volesti` relies on `Eigen` library [@eigen] for linear algebra but also support `MKL` optimizations [@mkl].
There are R and Python interfaces to `volesti` available [@Chalkis:2021, @Chalkis_dingo:2023]
although not all C++ features are exposed in those interfaces.

# Statement of need

High-dimensional sampling from multivariate distributions with Markov Chain Monte Carlo (MCMC)
algorithms is a fundamental problem with many applications in science and engineering [@Iyengar:1988;
@Somerville:1998; @Genz:2009; @Schellenberger:2009].
In particular, multivariate integration over a convex set and volume approximation of such sets
have accumulated a broad amount of effort over the last decades.
Nevertheless, those problems are computationally hard for general dimensions [@Dyer:1988].
MCMC algorithms have made remarkable progress efficiently solving the problems of sampling and
volume estimation of convex bodies while enjoying great theoretical guarantees [@Chen:2018; @Lee:2018;
@Mangoubi:2019]. However, theoretical algorithms cannot be applied efficiently as is to real-life
computations.
Therefore, practical algorithms have been designed by relaxing the theoretical guarantees and
applying new algorithmic and statistical techniques to perform efficiently while at the same time
meeting the requirements for high accuracy results [@Emiris:2014; @Cousins:2015; @Chalkis_volume:2023;
@Kook:2022].
The volume method of @Cousins:2015 and the sampling method of @Kook:2022 are available as `MATLAB`
packages.
All the methods mentioned above are implemented in `volesti` as well as the sampling algorithm by
@Chalkis_hmc:2023 and the algorithms for spectahedra by @Chalkis_spectra:2022.

To our knowledge `volesti` has been used in conducting research in electric power systems [@Venzke:2019]
and in probabilistic inference [Spallitta:2024].


# Acknowledgements

We would like to thank the contributors to the `volesti` library for their valuable contributions and
feedback.

# References