### import packages
import tskit
import tqdm
import math
import numpy as np
import multiprocessing as mp
import functools
import pickle
import os

### import C extension
import matrix

### epsilon function
def epsilon(x):
  N = x.shape[0]
  mean = x.mean()
  colmean = np.tile(x.mean(axis = 0), (N, 1))
  rowmean = colmean.T
  return x + mean - colmean - rowmean

### main function
def varGRM_C(trees, file = None, A = math.inf, B = 0, map_func = (lambda x:x), g = (lambda x: 1/(x*(1-x)))):
  if map_func == None:
    map_func = (lambda x:x)
  N = trees.num_samples
  egrm_C = matrix.new_matrix(N)
  egrm2_C = matrix.new_matrix(N)
  tmp1 = np.zeros(N)
  tmp2 = 0
  total_mu = 0
  pbar = tqdm.tqdm(total = trees.num_trees, 
                   bar_format = '{l_bar}{bar:30}{r_bar}{bar:-30b}',
                   miniters = trees.num_trees // 100,
                   file = file)
  
  for tree in trees.trees():
    l = map_func(tree.interval[1]) - map_func(tree.interval[0])
    if tree.total_branch_length == 0: continue
    for c in tree.nodes():
      descendants = list(tree.samples(c))
      n = len(descendants)
      if(n == 0 or n == N): continue
      t = max(0, min(A, tree.time(tree.parent(c))) - max(B, tree.time(c)))
      mu = l * t * 1e-8
      p = float(n/N)
      matrix.add_square(egrm_C, descendants, mu * g(p))
      matrix.add_square(egrm2_C, descendants, mu * np.power(g(p), 2) * np.power((1 - 2*p), 2))
      tmp1[descendants] += mu * np.power(g(p), 2) * (np.power(p, 2) - 2 * np.power(p, 3))
      tmp2 += mu * np.power(g(p), 2) * np.power(p, 4)
      total_mu += mu
    pbar.update(1)
  
  egrm = mat_C_to_array(egrm_C, N)
  egrm2 = mat_C_to_array(egrm2_C, N)
  
  egrm /= total_mu
  egrm2 /= total_mu
  tmp1 /= total_mu
  tmp2 /= total_mu
  
  e = np.reciprocal(np.random.poisson(lam=total_mu, size=10000).astype("float")).mean()
  egrm_final = epsilon(egrm)
  vargrm_final = e * (egrm2 + np.tile(tmp1, (N, 1)) + np.tile(tmp1, (N, 1)).transpose() + tmp2 - np.power(egrm_final, 2))
  
  pbar.close()
  return egrm_final, vargrm_final, total_mu