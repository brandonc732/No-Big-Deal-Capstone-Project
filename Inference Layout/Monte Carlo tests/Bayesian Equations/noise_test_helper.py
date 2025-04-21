import numpy as np
import pandas as pd

from scipy.stats import gamma


def run_noise_test(member):

    positions = member["positions"]
    confidences = member["confidences"]
    noise_mu = member["noise_mu"]
    noise_sigma = member["noise_sigma"]
    num_particles = member["num_particles"]
    prior_alpha = member["prior_alpha"]
    prior_beta = member["prior_beta"]
    radii = member["radii"]
    C = member["C"]
    
    # make noisy positions
    noisy_positions = positions + np.random.normal(loc=noise_mu, scale=noise_sigma, size=positions.shape)

    percentages = []
    variances = []

    # go through all particles
    # run through inference for each particle
    for i in range(num_particles):
        # get x and y positions
        x_pos = noisy_positions[:, i, 0]
        y_pos = noisy_positions[:, i, 1]

        # calculate movement by differences between positions
        x_dist = np.diff(x_pos)
        y_dist = np.diff(y_pos)

        # square all the distance movements
        x_dist_sqd = np.square(x_dist)
        y_dist_sqd = np.square(y_dist)

        # concatenate distances into one array
        dist_sqd_array = np.concatenate((x_dist_sqd, y_dist_sqd))
        
        # calculate posterior alpha and beta
        num_dists = len(dist_sqd_array)
        posterior_alpha = prior_alpha + (num_dists / 2)
        
        posterior_beta = (prior_beta + (sum(dist_sqd_array) / 2)) / C
        
        
        # get true radius
        true_radius = radii[i]

        # calculate cdp of true radius
        p_true_radius = gamma.cdf(true_radius, a=posterior_alpha, scale= 1 / posterior_beta)
        
        # append to percentages
        percentages.append(p_true_radius)

        # append to variances
        variances.append(gamma.var(a=posterior_alpha, scale= 1 / posterior_beta))

    percentages = np.array(percentages)

    results = {}
    for c in confidences:
        tail = (1 - c / 100) / 2
        results[c] = np.mean((percentages >= tail) & (percentages <= 1 - tail))
    
    results["mean_variance"] = np.mean(np.array(variances))
    
    return results


