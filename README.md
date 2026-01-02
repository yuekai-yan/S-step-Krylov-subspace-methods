# Fall 2025 Semester Project - Improving stability and performance of s-step Krylov subspace methods with randomization

## Overview
This project investigates deterministic and randomized orthogonalization schemes within the $s$-step GMRES framework. We systematically evaluate their impact on numerical accuracy, orthogonality preservation, and stability in finite-precision arithmetic, highlighting how randomization can improve robustness while maintaining computational efficiency.

## Folder Structure
```
├── +AOB                  # This directory contains inter-block orthogonalization methods
│   ├── CGS.m
│   ├── CGS2.m
│   ├── MGS.m
│   ├── RGS.m
│   ├── rCGS.m
│   ├── rCGS2.m
│   └── rMGS.m
├── +WB                   # This directory contains intra-block orthogonalization methods
│   ├── CGS.m
│   ├── CGS2.m
│   ├── Householder.m
│   ├── MGS.m
│   ├── RGS.m
│   ├── rCGS.m
│   ├── rCGS2.m
│   ├── rCholesky.m
│   └── rMGS.m
├── BGS_GMRES.m              # This script is used to implement BGS-GMRES
├── KrylovBasis              # This directory contains functions for generating monomial and Newton bases
│   ├── getRitzValues.m
│   └── mpk.m
├── OrthComp_BGS.m          # This script computes and records the relative residual and orthogonality loss  
                              in each iteration for different randomized orthogonalization schemes and 
                              step sizes s in BGS-GMRES 
├── OrthComparison.m        # This script computes and records the relative residual and orthogonality loss  
                              in each iteration for different randomized orthogonalization schemes and 
                              step sizes s in RBGS-GMRES 
├── RBGS_GMRES.m            # This script is used to implement RBGS-GMRES
├── SubspaceEmbedding       # This directory contains functions for generating subspace embeddings
│   ├── CountSketch.m
│   ├── Gaussian.m
│   ├── Rademacher.m
│   └── sparsesign.c
├── fig                       # This directory contains all the figures in the project
├── main.m                    # This script generates iteration-wise plots to compare deterministic and randomized
                                orthogonalization schemes of interest. It also supports studying the effect
                                of different step sizes s for a single or multiple deterministic/randomized schemes
├── plot4orthcomp.m           # This script plots heatmaps for the data obtained from `OrthComparison.m` and `OrthComp_BGS.m`, 
                                and allows selecting a suitable index and the corresponding data for visualization
├── plot_heatmap.m            # This script is used to plot the heatmap
└── test matrices             # This directory contains all the test matrices in the project
    ├── Si10H16.mat
    ├── SiH4.mat
    ├── SiO2.mat
    ├── getStartMatrix.m
    ├── laeuchli.m
    ├── mmread.m
    └── monomial_matrix.m
