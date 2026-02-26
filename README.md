# MemRBC

MemRBC is a result of open science research. It demonstrates the use of hi-level apps in a light-weight object oriented approach in R, keeping metadata and lineage recorded systematically for reproducibility and self-documentation. MemRBC is the R package for modeling red blood cell shape. It is based on a spectral model, using spherical harmonics functions for a continuous parametrerization of 3D shape (of genus 0). Besides usual Helfrich, area difference elasticity energy and spontaneous curvature, the model includes a shear-elastic network component to model spectrin network mechanics in the same numerical scheme.

MemRBC gives a set of high-level routines to the user - the Apps - to perform typical modeling tasks like energy minimization or Montecarlo dynamics simulation. From the codes behind these Apps, the user may develop his/her own extensions and modifications, and take insight in the backbone code infrastructure.

The philosophy of MemRBC data model is to keep as much information as possible on the performed steps within the MemRBC object. Besides a App-command history, energy values, iteration counts and more, total processing time is recorded. Since MemRBC is a class, it may use plot, print, update and further existing generics.

MemRBC makes use of Rcpp, for faster vectorizing numerics in the energy, gradient and Hessian computations. However, the original R implementations are also available. The C++-code is shipped within the MembraneRBC.R file for inline compilation. One may easily modify and make experiments with parallelisation, e.g. using more OpenMP.

The idea of continuous parametrization originates from the dissertation thesis\
Frickenhaus, S. (1999) Modelling of Lipid Translocation - Application to Red Blood Cell Membrane. Shaker Verlag , Aachen, 1999, ISBN: 978-3-8265-6366-9.

### Updates
Updates are found on Github; latest release is 0.4.0

### Acknowledgement

Thanks go to Khaled Khairy for pointing at the spherical harmonics utilities and a collection of shape coefficients, and Wolfgang Bosch for comments on the original spherical harmonics code. Mathis Pink is acknowledged for discussions on soft constraints, genetic algorithms and alternatives to MMC, like NUTS.

### Licensing

The code is licensed under GPL 3.
