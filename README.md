# No Big Deal Capstone Project
Project Title: **Improved Nanoparticle Size Posterior Estimation by Markov Chain Monte Carlo Inference with Location Error Characterization**

Team No Big Deal:
- Brandon Chestnut
- Soham Kapadia
- Aidan Shea
- Harshit Singh

**Soham's Simulation:** [GitHub](https://github.com/soso1407/nanoparticle-tracking-simulation)


# Getting Started

- Since we are not allowed to upload NIST images to the public GitHub, you need to copy the folder of NIST test data supplied by Dr Pintar into the `test_data` folder.

- We have provided the Anaconda Python environment used to run the project in the `anaconda environment` folder. All the libraries are pretty standard except for a few of the backend accelerations used for PyMC, which could potentially lead to issues across platforms.
   - PyMC also offers CUDA GPU acceleration, but only on Linux, so it was not used here.

Otherwise, all Python components in the respository should be plug and play. The beads program could pose some trouble, but more detail for this is provided in [beads readme](https://github.com/brandonc732/No-Big-Deal-Capstone-Project/blob/main/beads/readme.md)

# Main notebooks

This folder includes the primary large notebooks of the project. Each Jupyter notebook should outline its process with markdown sections and comments.

**note:** it's highly recommended to collapse jupyter notebook cell sections for readability 

## Simulation.ipynb

author : Brandon

With the project goal of modeling or characterizing the location error within the NIST testing data, `Simulation.ipynb` Creates a simulation fitted to mimic the appearance of the test images. Using strong TrackPy identifications in the test data, the simulation:

- Mimics background noise by fitting a guassian distribution to test image background pixels (0 - 0.9999 percentile)
- Mimics particle appearance by fitting 2D guassian curves to each TrackPy identification. The resulting x_sigmas and y_sigmas are sampled according to the signal of their identification.
- Mimics particle brightness change by having simulated particle brightnesses use a random walk similar to that seen in the TrackPy identifications

The steps to this process along with other image details are outlined with markdown sections and plots throughout the notebook.

Overall, the simulation was able to produce images that look similar to the test data both in particle shape and sporadic particle behavior. However, just recreating the fitted 2D guassian curves in the simulation images did not induce any significant errors in subsequent TrackPy identifications. This is detailed in final section of the notebook.

Further changes would need to be made to qualify the aspects of the testing data that induce localization error such as pixel deviation from the Gaussian blur. This is implemented in `SimulationUpdate.ipynb`.

Example [simulation gif](readme_images/sim_gamma.gif) (not allowed to have NIST video in github to compare).


## SimulationUpdate.ipynb

Soham and Aidan's changes to improve original simulation notebook.

## Knob Framework & Grid‑Search Calibration (Soham’s Contribution)

- Introduced three global tuning parameters to bridge simulation and real data:
  1. **noise_amp**: scales the background‑noise grain  
  2. **sigma_mul**: rescales the per‑particle PSF width (bright spots sharpened when < 1, dim ones blurred when > 1)  
  3. **grad_x, grad_y**: linear illumination gradient across the frame

- Performed a discrete grid search over all 81 combinations of  
  noise_amp ∈ {0.5, 1.0, 1.5}, sigma_mul ∈ {0.8, 1.0, 1.2}, (grad_x, grad_y) ∈ {−0.01, 0, 0.01}²  
  on the first ten frames (noise‑free branch) to find the triple minimizing pixel‑wise RMS error.

- Locked in the best parameters (noise_amp=1.0, sigma_mul=0.8, background_grad=(+0.01, −0.01)) for all subsequent high‑fidelity renders and MCMC inference.


Aidan created a function called "fit_pixel_noise_kde" that builds a realistic pixel-level noise model by extracting residuals from real particle image patches. For each patch it fits a 2d gaussian, and computes the residual noise. These residuals are then used to fit a kernel density estimate (KDE) to find the empirical distribution of pixel noise. The kde is then used in the simulation to add any imperfections of camera based noise.

## CheatMethod.ipynb

author: Harshit

Built off of MCMC code in `Inference Layout.ipynb`
This Jupyter notebook provides a calibration framework for nanoparticle tracking by implementing a "Cheat Method" error modeling approach. The core functionality simulates Gaussian localization noise on known (true) particle trajectories, systematically adjusting the injected noise magnitude until the mean estimated particle radius, obtained via the Stokes-Einstein relation, matches the nominal radius. The workflow includes:

Synthetic trajectory generation with precise control over noise characteristics,

Localization error injection via Gaussian perturbations,

Visualization of particle paths and uncertainty rings using TrackPy and Matplotlib,

Manual tuning of pixel-level noise parameters based on ground-truth recovery.
This notebook is designed to benchmark and recalibrate particle tracking pipelines where ground truth is available, offering a robust tool to fine-tune localization uncertainty estimates and ensure physically meaningful radius recovery.


# Folders

Each folder includes its own readme describing its contents and purpose.

## Inference Layout

author : Brandon

Contains all work for deriving inference methods through ground truth demonstrations. This includes `Inference Layout.ipynb` where we create initial working versions of Bayesian Equation and MCMC posterior inference methods along with Monte Carlo demonstration folders where we test their long term performance

## beads

Contains all work done with beads.exe from [![DOI](https://img.shields.io/badge/DOI-10.1214%2F09--AOAS299-blue)](https://doi.org/10.1214/09-AOAS299)

Ultimately couldn't use it for the project. More detail in [beads readme](https://github.com/brandonc732/No-Big-Deal-Capstone-Project/blob/main/beads/readme.md)

## test data

location to insert test data from folder Dr. Pintar shared. Should include `Fluorescent_99_nm_polystyrene_in_saline` folder and `particle_tracking_video.avi` video

## anaconda environment

folder with yaml file for Anaconda environment to run project.
