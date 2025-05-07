Author : Brandon

# Overview

This folder contains the two folders used for demonstrating expected performance the two solution approaches by running inference for a large number of particles (Monte Carlo setup).

Each notebook is sectioned off into different Monte Carlo tests with simple descriptions of each. For detail on how the method works or is derived, see [Inference Layout.ipynb](https://github.com/brandonc732/No-Big-Deal-Capstone-Project/blob/main/Inference%20Layout/Inference%20Layout.ipynb)


# Multiprocessing and RAM

To help with computing time for large tests such as running Monte Carlo for different location error variances, both notebooks use `Pool` from Python's `multiprocessing` module. This is a very basic version of multiprocessing that spawns a new Python process for each instance in the population.

Since my desktop has 64 GB of RAM, I didn't bother optimizing the multiprocessing method. So when running intensive tasks such as the MCMC tests, I get RAM useages up to 40 gigabytes.

**note:** If RAM usage prevents the code from running, specifically with MCMC, change variables such:
- num_chains
- num_particles
- num_sigmas
- num_cores
- change draws or tune within [no_noise_MCMC_helper.py](MCMC/no_noise_MCMC_helper.py) or [normal_noise_MCMC_helper.py](MCMC/normal_noise_MCMC_helper.py)
  - Will be the most effective in reducing RAM usage and computation, but of course results in less posterior exploration.



# Bayesian Equations folder

Includes Monte Carlo demonstrations for posterior estimation with Gamma equations described in [Inference Layout.ipynb](https://github.com/brandonc732/No-Big-Deal-Capstone-Project/blob/main/Inference%20Layout/Inference%20Layout.ipynb) in the following situations:
- direct ground truth locations
- direct ground truth locations with set guassian noise
- credible set coverage of previous Monte Carlo test along a range of gaussian noise standard deviation 

**note:** since inference is straightforward equations, Python can easy handle over 100,000 particle counts


# MCMC folder

This folder demonstrates similar scenarios for MCMC, however, we focus on whether MCMC is able to account for location error variance in our estimation, not just how uninformed MCMC estimation is affected.

Includes Monte Carlo demonstrations for posterior estimation using MCMC as described in [Inference Layout.ipynb](https://github.com/brandonc732/No-Big-Deal-Capstone-Project/blob/main/Inference%20Layout/Inference%20Layout.ipynb) in the following situations:
- direct ground truth locations
- direct ground truth locations with set guassian noise
  - perfect error variance knowledge
- credible set coverage of previous Monte Carlo test along a range of gaussian noise standard deviation
  - perfect error variance knowledge
- credible set coverage of previous Monte Carlo test along a range of gaussian noise standard deviation
  - 10% error variance underestimation

**note:** since MCMC inference is  computationally expensive, I only tested particle counts of around 100 to keep run times reasonable


## Numpyro

The standard NUTS sampler in PyMC is simply too computationally expensive to use in Monte Carlo demonstrations. 

So I employed the `numpyro` NUTS sampler, which is significantly faster. However, it is known to have random compatability issues. For example, it wasn't compatible with some aspect of widgets in my version of Jupyter notebook and simply wouldn't run. 

The [Anaconda environment](https://github.com/brandonc732/No-Big-Deal-Capstone-Project/blob/main/anaconda%20environment/Capstone_environment.yaml) within this github should have all compatability issues sorted, but if you run into to anything, try changing the code to use the standard sampler instead to rule out numpyro issues.

**note:** Numpyro offers CUDA GPU acceleration for Linux















