### importing packages
#import sys
#import time
#import datetime
import msprime
#import tskit
import math
import numpy as np
#import pandas as pd
#import os

### parameters
out = "tmp"
log_file = "tmp.log"

l = 32000000
N = 1000
mutation_rate = 1e-8
recomb_rate = 1e-8

cas_ratio = 0.1
obs_ratio = 0.2
Alpha = -0.5
Beta = 1
h2g = 1.0

### out-of-Africa parameters
N_B = 2100
N_EU0 = 1000
N_AF = 12300
N_A = 7300
r_EU = 0.004
generation_time = 25
T_EU_AS = 21.2e3 / generation_time
T_B = 140e3 / generation_time
T_AF = 220e3 / generation_time
N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)

### parameters
def simulate_hapdata(l = l, N = N, mutation_rate = mutation_rate, recomb_rate = recomb_rate):
    population_configurations = [msprime.PopulationConfiguration(sample_size = N, initial_size = N_EU, growth_rate = r_EU)]
    
    demo_events = [msprime.PopulationParametersChange(time=T_EU_AS, initial_size = N_B, growth_rate=0),
                   msprime.PopulationParametersChange(time=T_B, initial_size = N_AF, growth_rate=0),
                   msprime.PopulationParametersChange(time=T_AF, initial_size = N_A, growth_rate=0)]
    
    tree_sequence = msprime.simulate(length = l, population_configurations = population_configurations, 
                                     recombination_rate = recomb_rate, mutation_rate = mutation_rate, 
                                     demographic_events = demo_events)
    
    genotype_matrix = tree_sequence.genotype_matrix().astype("uint16")
    MAFs = genotype_matrix.mean(axis = 1); MACs = np.array(list(map(min, MAFs, 1-MAFs)))
    variants = np.array([v.position for v in tree_sequence.variants()])
    M = len(variants)
    
    return {"tree_sequence":tree_sequence, "genotype_matrix":genotype_matrix, 
            "MAFs":MAFs, "MACs":MACs, "variants":variants, "M":M}

def simulate_phenotypes(hapdata, h2g = h2g, cas_ratio = cas_ratio, Alpha = Alpha):
    genotype_matrix = hapdata["genotype_matrix"]
    N = genotype_matrix.shape[1]
    M = hapdata["M"]
    MAFs = hapdata["MAFs"]
    
    M_cas = int(round(M * cas_ratio))
    weights = np.multiply(np.power(MAFs, Alpha), np.power(1-MAFs, Alpha))
    weights /= weights.sum()
    cass = np.random.choice(range(M), M_cas, replace = False, p = weights); cass.sort()
    
    X_cas = np.transpose(genotype_matrix[cass])
    Z_cas = X_cas - np.mean(X_cas, axis=0)
    Z_cas = Z_cas / np.std(Z_cas, axis=0)
    
    betas = np.random.normal(scale = np.sqrt(h2g / M_cas), size = (M_cas, 1))
    
    G = np.dot(X_cas, betas); G.shape = (N, )
    s2g = np.var(G, ddof = 1)
    s2e = s2g * (1/h2g - 1)
    
    e = np.random.normal(scale=np.sqrt(s2e), size=N)
    y = G + e
    
    maternals = np.arange(0, N-1, 2)
    paternals = np.arange(1, N, 2)
    y_diploid = y[maternals] + y[paternals]
    
    return {"M_cas":M_cas, "cass":cass, "Z_cas":Z_cas, "betas":betas, "y":y, "y_diploid":y_diploid}

def simulate_observations(hapdata, obs_ratio = obs_ratio, Beta = Beta):
    genotype_matrix = hapdata["genotype_matrix"]
    variants = hapdata["variants"]
    MAFs = hapdata["MAFs"]
    MACs = hapdata["MACs"]
    
    idx_5 = MACs >= 0.005
    M_5 = idx_5.sum()
    MAFs_5 = MAFs[idx_5]
    
    variants_ceil = np.ceil(variants)
    overlapped = np.insert(variants_ceil[:-1] == variants_ceil[1:], 1, False)
    non_overlapped = np.logical_not(overlapped)
    
    observable =  np.logical_and(idx_5, non_overlapped)
    M_observable = observable.sum()
    
    weights = np.multiply(np.power(MAFs, Beta), np.power(1-MAFs, Beta))
    weights = np.multiply(weights, observable)
    weights /= weights.sum()
    
    M_obs = int(round(M_observable * obs_ratio))
    obss = np.random.choice(np.where(observable)[0], M_obs, replace = False, p = weights[observable]); obss.sort()
    
    X_obs = np.transpose(genotype_matrix[obss])
    Z_obs = X_obs - np.mean(X_obs, axis=0)
    Z_obs = Z_obs / np.std(Z_obs, axis=0)
    
    return {"M_5":M_5, "M_obs":M_obs, "obss":obss, "Z_obs":Z_obs}

def simulate(l = l, N = N, mutation_rate = mutation_rate, recomb_rate = recomb_rate,
            h2g = h2g, cas_ratio = cas_ratio, Alpha = Alpha,
            obs_ratio = obs_ratio, Beta = Beta):
    hapdata = simulate_hapdata()
    phenotypes = simulate_phenotypes(hapdata)
    observations = simulate_observations(hapdata)
    
    return {"hapdata":hapdata, "phenotypes":phenotypes, "observations":observations}



    

    
