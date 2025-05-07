Author : Brandon

# Overview

This folder contains all the Jupyter notebooks used in deriving solution approaches with Bayesian Equations and Markov Chain Monte Carlo (MCMC).

In `Inference Layout.ipynb`, I detail the math and distributions within the markdown cells to explain how posteriors are derived.

`Inference Layout.ipynb` includes:
- Brownian motion simulation equations
- Bayesian posterior distribution through gamma distribution
   - 2d and 3d case, though 3d is not pursued in the project
- MCMC with PyMC (no location error)
- MCMC with PyMC (undefined location error variance)
  - demonstrates that MCMC fails with undefined location error variance
- MCMC with PyMC on 20 particle simulation
  - demonstrates that MCMC works when we can define the location error variance. In this case, it's defined by comparing TrackPy locations against the [ground truth csv](https://github.com/brandonc732/No-Big-Deal-Capstone-Project/blob/main/Inference%20Layout/20_particles/simulation_v3_tracks.csv).

We investigate how error affects the credible set coverage of each method in the `Monte Carlo tests` folders







