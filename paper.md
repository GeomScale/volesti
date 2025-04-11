---
title: 'volesti: A C++ library for sampling and volume computation on convex bodies'
tags:
  - C++
  - geometry
  - randomization
  - Monte Carlo methods
  - convexity
authors:
  - name: Apostolos Chalkis
    orcid: 0000-0002-4628-1907
    equal-contrib: true
    affiliation: "2, 4" # (Multiple affiliations must be quoted)
  - name: Vissarion Fisikopoulos
    orcid: 0000-0002-0780-666X
    corresponding: true # (This is how to denote the corresponding author)
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: "1, 4"
  - name: Marios Papachristou
    orcid: 0000-0002-1728-0729
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 5
  - name: Elias Tsigaridas
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: "3, 4"
affiliations:
 - name: National & Kapodistrian University of Athens, Greece
   index: 1
 - name: Quantagonia, Germany
   index: 2
 - name: Inria Paris and IMJ-PRG, Sorbonne Université, France
   index: 3
 - name: GeomScale, Greece
   index: 4
 - name: Cornell University, U.S.A.
   index: 5
date: 11 March 2024
bibliography: paper.bib

---

# Summary

Sampling from (constrained) high-dimensional distributions and volume approximation of convex
bodies are fundamental operations that appear in optimization, finance,
engineering, artificial intelligence, and machine learning.
We present `volesti`, a C++ library that delivers efficient implementations of  state-of-the-art, mainly randomized, algorithms
to sample from general logarithmically concave (or log-concave) distributions.
Based on these routines, we can estimate the volume of convex bodies in high dimensions,
round them, and compute multidimensional integrals over them.
The backbone of our library consists of Monte Carlo algorithms,
which are randomized algorithms, the output of which can be incorrect with (usually very small) error probability; thus, we also provide several
high-dimensional statistical tests to certify and verify the output.

The focus of `volesti` is scalability in high dimensions,
that, depending on the problem at hand, could range from hundreds to thousands of dimensions.
Another novelty is the ability to handle a variety of different inputs
for the constrained support of the various distributions.
`volesti` supports three different types of polyhedra [@Ziegler:1995], spectrahedra [@Ramana:1999],
and general non-linear convex objects.

`volesti` relies on `Eigen` library [@eigen] for linear algebra but also supports `MKL` optimizations [@mkl].
There are R [@Chalkis:2021] and Python [@Chalkis_dingo:2023] interfaces available.

# Statement of need

High-dimensional sampling from multivariate distributions with Markov Chain Monte Carlo (MCMC)
algorithms is a fundamental problem with many applications in science and engineering [@Iyengar:1988; @Somerville:1998; @Genz:2009; @Schellenberger:2009].
In particular, multivariate integration over a convex set as well as the volume approximation of convex sets have garnered significant attention from theorists and engineers  over the last decades.
Nevertheless, these problems are computationally hard for general dimensions [@Dyer:1988].
MCMC algorithms made remarkable progress and their use allowed us to efficiently tackle the problems of sampling and volume estimation of convex bodies in theory,
by the introduction of (rigorous) theoretical guarantees [@Chen:2018; @Lee:2018;
@Mangoubi:2019].
Unfortunately, these theoretical guarantees of the MCMC algorithms
do not extend in a straightforward manner to efficient implementations able to attack problems coming from real-life computations.
Therefore, we witnessed the birth of efficient in practice MCMC algorithm
that relax the theoretical guarantees and employ new algorithmic and statistical techniques
to be amenable to efficient implementations.
Remarkably, these algorithms, and the corresponding implementations,
also meet the requirements for high accuracy results
[@Emiris:2014; @Cousins:2015; @Chalkis_volume:2023; @Kook:2022];
however several existing published methods are only available as part of propertiary packages (MATLAB) [@Cousins:2015; @Kook:2022].

Our open-source package -- `volesti` -- offers all of the aforementioned functionality, together with the support of sampling from general log-concave densities [@Chalkis_hmc:2023], and uniform sampling from spectrahedra [@Chalkis_spectra:2022].

Our implementation:

1. supports various sampling techniques based on geometric walks; roughly speaking these are a continuous version of MCMC algorithms, such as Billard walk, Hamiltonian walk and others,
2. gives the user the ability to sample from  various distributions, like uniform, exponential, Gaussian, and general log-concave densities,
3. allows to consider the distributions constrained in various convex domains, such as hypercubes, zonotopes, general polytopes (defined either as a set of linear inequalities or as a convex hull of a pointset), spectrahedra (feasible sets of semidefinite programs), and
4. can perform volume computations, integration, and solve problems from real life applications in very high dimensions.

# Impact

`volesti` has been used extensively in various research and engineering projects coauthored by the authors of this paper.
In particular, for the problem of sampling the flux space of metabolic networks
we were able to sample from the  most complicated human metabolic network accessible today, Recon3D [@cftz-socg021],
to model financial crises [@ccef-crises-j],
to detect low volatility anomalies in stock markets [@bcft-aistats-23],
   to introduce randomized control in asset pricing and portfolio performance evaluation [@bcft-arxiv-24]), and also to sample from (and compute the volume of) spectrahedra [@Chalkis_spectra:2022], the feasible regions of semidefinite programs.

`volesti` has also been used by other research teams in conducting research in electric power systems [@Venzke:2019], for problems in probabilistic inference [@Spallitta:2024],
to perform resource analysis on programs [@pham-phd-2024];
and also for more theoretical and mathematical challenges, like the computation of topological invariants [@co-alenex-2021] and persistent homology [@vm-fods-2022].

# Acknowledgements

We would like to thank the contributors to the `volesti` library for their valuable contributions and feedback.
MP was partially supported by a Cornell University Fellowship, a grant from the A.G. Leventis Foundation, a grant from the Gerondelis Foundation, and a LinkedIn Ph.D. Fellowship.

# References
