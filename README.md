# No Big Deal Capstone Project
Project Title: **Improved Nanoparticle Size Posterior Estimation by Markov Chain Monte Carlo Inference with Location Error Characterization**

Team No Big Deal:
- Brandon Chestnut
- Soham Kapadia
- Aidan Shea
- Harshit Singh

## Project overview

asdf



# Getting Started

- Since we are not allowed to upload NIST images to the public GitHub, you need to copy the folder of NIST test data supplied by Dr Pintar into the `test_data` folder.

- We have provided the anaconda Python environment used to run the project in the `anaconda environment` folder. All the libraries are pretty standard except for a few of the backend accelerations used for PyMC, which could potentially lead to issues across platforms.
   - PyMC also offers CUDA GPU acceleration, but only on Linux, so it was not used here.

Otherwise everything Python related in the respository should be plug and play. The beads program could pose some trouble, but more detail for this is provided in [beads readme](https://github.com/brandonc732/No-Big-Deal-Captsone-Project/blob/main/beads/readme.md)

# Main notebooks

This folder includes the primary large notebooks of the project. Each Jupyter notebook should outline its process with markdown sections and comments.

**note:** it's highly recommended to collapse jupyter notebook cell sections for readability 

## Simulation.ipynb

author : Brandon

With the project goal of modeling or characterizing the location error within the NIST testing data, `Simulation.ipynb` Creates a simulation fitted to mimick the appearance of the test images. Using strong TrackPy identifications in the test data, the simulation:

- Mimicks background noise by fitting a guassian distribution to test image background pixels (0 - 0.9999 percentile)
- Mimicks particle appearance by fitting 2D guassian curves to each TrackPy identification. The resulting x_sigmas and y_sigmas are sampled according to the signal of their identification.
- Mimicks particle brightness change by having simulated particle brightnesses use a random walk similar to that seen in the TrackPy identifications

The steps to this process along with other image details are outlined with markdown sections and plots throughout the notebook.

Overall, the simulation was able to produce images that look similar to the test data both in particle shape and sporatic particle behaivor. However, just recreating the fitted 2D guassian curves in the simulation images did not induce any significant errors in subsequent TrackPy identifications. This is detailed in final section of the notebook.

Further changes would need to be made to qualify the aspects of the testing data that induce localization error such as pixel deviation from the gaussian blur. This is implemented in `SimulationUpdate.ipynb`.

Example [simulation gif](readme_images/sim_gamma.gif) (not allowed to have NIST video in github to compare).


## SimulationUpdate.ipynb

Soham and Aidan's changes to improve original simulation notebook.

Describe changes made and how to run it since it's setup for google colab

## CheatMethod.ipynb

author: Harshit

Built off of MCMC code in `Inference Layout.ipynb` to ....


# Folders

Each folder includes its own readme describing its contents and purpose.

## Inference Layout

author : Brandon

Contains all work with deriving inference methods through ground truth demonstrations. This includes `Inference Layout.ipynb` where we create initial working versions of Bayesian Equation and MCMC posterior infernence methods along with Monte Carlo demonstration folders where we test their long term performance

## beads

Contains all work done with beads.exe from [![DOI](https://img.shields.io/badge/DOI-10.1214%2F09--AOAS299-blue)](https://doi.org/10.1214/09-AOAS299)

Ultimately couldn't use it for the project. More detail in [beads readme](https://github.com/brandonc732/No-Big-Deal-Captsone-Project/blob/main/beads/readme.md)

## test data

location to insert test data from folder Dr. Pintar shared. Should include `Fluorescent_99_nm_polystyrene_in_saline` folder and `particle_tracking_video.avi` video

## anaconda environment

folder with yaml file for anaconda environment to run project.
