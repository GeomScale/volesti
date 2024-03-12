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
 - name: Inria Paris and IMJ-PRG, Sorbonne Universit\`e
   index: 3
 - name: GeomScale
   index: 4
 - name: Cornell University
   index: 5
date: 11 March 2024
bibliography: paper.bib

---

# Summary

Sampling from (constrained) high-dimensional distributions and volume approximation of convex
bodies are fundamental operations that appear in optimization, finance,
engineering, artificial intelligence, and machine learning.
We present `volesti`, a C++ library that delivers efficient implementations of  state-of-the-art, mainly randomized, algorithms
to sample from general logconcave distributions.
Based on these routines can estimate the volume of convex bodies in high dimensions, 
round them and also compute multidimensional integrals over them.
The backbone of our library consists of Monte-Carlo algorithms,
that are randomized algorithms, the output of which can be incorrect with (usually very small) error probability; thus, we also provide several 
high-dimensional statistical tests to certify and verify the output. 

The focus of `volesti`' is scalability in high dimensions, 
that, depending on the problem at hand, could be in the order of hundreds or thousands dimension.
Another novelty is the ability to handle a variety of different inputs
for the constrained support of the various distributions.
`volesti` supports three different types of polyhedra [@Ziegler:1995], spectrahedra [@Ramana:1999]
and general non-linear convex objects.

`volesti` relies on `Eigen` library [@eigen] for linear algebra but also support `MKL` optimizations [@mkl].
There are R [@Chalkis:2021] and Python [@Chalkis_dingo:2023] interfaces available;
alas not all C++ functionality is available in through these interfaces.

# Statement of need

High-dimensional sampling from multivariate distributions with Markov Chain Monte Carlo (MCMC)
algorithms is a fundamental problem with many applications in the whole spectrum of  science and engineering [@Iyengar:1988;
@Somerville:1998; @Genz:2009; @Schellenberger:2009].
In particular, multivariate integration over a convex set
as well as the volume approximation of convex sets
have accumulated a huge amount of effort from theorists and engineers  over the last decades.
Nevertheless, these problems are computationally hard for general dimensions [@Dyer:1988].
MCMC algorithms made remarkable progress 
and their use allowed us to efficiently tackle the problems of sampling and
volume estimation of convex bodies in theory,
by the introduction of (ragher sharp) theoretical guarantees [@Chen:2018; @Lee:2018;
@Mangoubi:2019].
Unfortunately, these theoretical guarantees of the MCMC algorithms 
do not extend in an straightforward manner to efficient implementations able to attack problems coming from real-life computations.
Therefore, we witnessed the birth of efficient in practice MCMC algorithm
that they relax the theoretical guarantees and
and employ new algorithmic and statistical techniques 
to be amenable to efficient implementations.
Remarkably, these algorithms, and the corresponding implementations,
also meet the requirements for high accuracy results 
[@Emiris:2014; @Cousins:2015; @Chalkis_volume:2023; @Kook:2022].
Let us mention that the volume algorithm of @Cousins:2015 and the sampling method of @Kook:2022 are available as `MATLAB`
packages.

All aforementioned algorithms and techniques are available in `volesti`
along with the  sampling algorithm by
@Chalkis_hmc:2023 and the algorithms for spectrahedra by @Chalkis_spectra:2022.

The efficient implementations of `volesti  
(i) suport various sampling techniques based on geometric walks, roughly speaking these are a continuous version of MCMC algorithms, like Billard walk, Hamiltonian walk and other,
(ii) give us the ability to sample from  various distributions, like uniform, log-concave, exponential, and Gaussian,
(iii) allows to consider the distributions 
    constrained in various convex domains, like hypercubes, zonotopes, general polytopes (in H and V representations), spectrahedra,
    and (iv) can perform volume computations, integration, and solve problem from real life applications in very high dimensions.

    
    
We use `volesti` extensively in various research and engineering directions that we pursue.
In particular, for the problem of sampling the flux space of metabolic networks
we were able to sample from the  most complicated human metabolic network accessible today, Recon3D [@cftz-socg021],
we use to model financial crises [@ccef-crises-j],
to detect low volatility anomalies in stock markets [@bcft-aistats-23],
   to introduce randomized control in asset pricing and portfolio performance evaluation [@bcft-arxiv-24]), but also to sample from (and compute the volume of) spectrahedra [@Chalkis_spectra:2022], the feasible regions of semidefinite programs. 
    
Even more, `volesti` has been used in conducting research in electric power systems [@Venzke:2019], for problems 
in probabilistic inference [@Spallitta:2024],
to perform resource analysis on programs [@pham-phd-2024]; 
but also to more theoretical and mathematical challenges, like the computation of topological invariants [@co-alenex-2021]
    and persistent homology [@vm-fods-2022].

# Acknowledgements

We would like to thank the contributors to the `volesti` library for their valuable contributions and
feedback.

# References
