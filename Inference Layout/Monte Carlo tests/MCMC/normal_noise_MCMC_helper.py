import numpy as np

import pymc as pm
import arviz as az


def run_particle_range(member):
    particle_nums = member["particle_nums"] 
    positions = member["positions"]
    radii = member["radii"]
    C = member["C"]
    num_chains = member["num_chains"]
    true_noise_mean = member["true_noise_mean"]
    true_noise_sigma = member["true_noise_sigma"]
    estimated_noise_sigma = member["estimated_noise_sigma"]
    
    # make noisy positions
    noisy_positions = positions + np.random.normal(loc=true_noise_mean, scale=true_noise_sigma, size=positions.shape)

    percentages = []
    variances = []

    # go through particles in assigned range
    for i in particle_nums:
        # get x and y positions
        x_pos = noisy_positions[:, i, 0]
        y_pos = noisy_positions[:, i, 1]

        # calculate movement by differences between positions
        x_dist = np.diff(x_pos)
        y_dist = np.diff(y_pos)

        dist_arr = np.concatenate((x_dist, y_dist))    
    
        # get true radius
        true_radius = radii[i]

        sigma_e2 = estimated_noise_sigma**2

        # Define model
        with pm.Model() as model:
            #tau = pm.Gamma("tau", alpha=prior_alpha, beta=prior_beta)
            # Prior for process variance sigma^2 (can adjust based on domain knowledge)
            sigma2 = pm.HalfCauchy("sigma2", beta=2)
            
            tau = pm.Deterministic("tau", 1 / sigma2)
            
            # Total standard deviation of observed deltas
            total_std = pm.Deterministic("total_std", pm.math.sqrt(sigma2 + 2 * sigma_e2))

            # Likelihood
            obs = pm.Normal("obs", mu=0, sigma=total_std, observed=dist_arr)

            # Use NUTS sampler
            trace = pm.sample(draws=20000, tune=10000, chains=num_chains, target_accept=0.9, return_inferencedata=True, nuts_sampler="numpyro", progressbar=False)

        # extract radius samples
        tau_samples = trace.posterior["tau"].values.flatten()
        rad_samples = tau_samples * C
        
        # calculate cdf of true radius
        p_true_radius = np.mean(rad_samples < true_radius)

        # append to percentages
        percentages.append(p_true_radius)

        # append to variances
        variances.append(np.var(rad_samples))
    
    percentages = np.array(percentages)
    variances = np.array(variances)

    return [percentages, variances]






def run_variance_option(member):
    num_particles = member["num_particles"] 
    positions = member["positions"]
    radii = member["radii"]
    C = member["C"]
    num_chains = member["num_chains"]
    true_noise_mean = member["true_noise_mean"]
    true_noise_sigma = member["true_noise_sigma"]
    estimated_noise_sigma = member["estimated_noise_sigma"]
    confidences = member["confidences"]
    
    # make noisy positions
    noisy_positions = positions + np.random.normal(loc=true_noise_mean, scale=true_noise_sigma, size=positions.shape)

    percentages = []
    variances = []

    # go through particles in assigned range
    for i in range(num_particles):
        # get x and y positions
        x_pos = noisy_positions[:, i, 0]
        y_pos = noisy_positions[:, i, 1]

        # calculate movement by differences between positions
        x_dist = np.diff(x_pos)
        y_dist = np.diff(y_pos)

        dist_arr = np.concatenate((x_dist, y_dist))    
    
        # get true radius
        true_radius = radii[i]

        sigma_e2 = estimated_noise_sigma**2

        # Define model
        with pm.Model() as model:
            #tau = pm.Gamma("tau", alpha=prior_alpha, beta=prior_beta)
            # Prior for process variance sigma^2 (can adjust based on domain knowledge)
            sigma2 = pm.HalfCauchy("sigma2", beta=2)
            
            tau = pm.Deterministic("tau", 1 / sigma2)
            
            # Total standard deviation of observed deltas
            total_std = pm.Deterministic("total_std", pm.math.sqrt(sigma2 + 2 * sigma_e2))

            # Likelihood
            obs = pm.Normal("obs", mu=0, sigma=total_std, observed=dist_arr)

            # Use NUTS sampler
            trace = pm.sample(draws=20000, tune=10000, chains=num_chains, target_accept=0.9, return_inferencedata=True, nuts_sampler="numpyro", progressbar=False)

        # extract radius samples
        tau_samples = trace.posterior["tau"].values.flatten()
        rad_samples = tau_samples * C
        
        # calculate cdf of true radius
        p_true_radius = np.mean(rad_samples < true_radius)

        # append to percentages
        percentages.append(p_true_radius)

        # append to variances
        variances.append(np.var(rad_samples))


    percentages = np.array(percentages)

    confidences = member["confidences"]

    results = {}
    for c in confidences:
        tail = (1 - c / 100) / 2
        results[c] = np.mean((percentages >= tail) & (percentages <= 1 - tail))
    
    results["mean_variance"] = np.mean(np.array(variances))
    
    return results
