import numpy as np

import pymc as pm
import arviz as az



def run_no_noise_MCMC(member):
    particle_nums = member["particle_nums"] 
    positions = member["positions"]
    radii = member["radii"]
    C = member["C"]
    num_chains = member["num_chains"]

    percentages = []
    variances = []


    for i in particle_nums:
        # get x and y positions
        x_pos = positions[:, i, 0]
        y_pos = positions[:, i, 1]

        # calculate movement by differences between positions
        x_dist = np.diff(x_pos)
        y_dist = np.diff(y_pos)

        dist_arr = np.concatenate((x_dist, y_dist))    
    
        # get true radius
        true_radius = radii[i]


        # Define model
        with pm.Model() as model:
            #tau = pm.Gamma("tau", alpha=prior_alpha, beta=prior_beta)
            #sigma = pm.Deterministic("sigma", 1 / pm.math.sqrt(tau))
        
            # suggested method for non informative prior (works way better)
            sigma = pm.HalfCauchy("sigma", beta=5)
            tau = pm.Deterministic("tau", 1 / sigma**2)
            
            obs = pm.Normal("obs", mu=0, sigma=sigma, observed=dist_arr)
            
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
    

























