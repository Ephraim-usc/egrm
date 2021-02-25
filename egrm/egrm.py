### import packages
import tskit
import tqdm
import math
import numpy as np


### import C extension
import matrix


def mat_C_to_array(mat_C, N):
  buffer = matrix.export_matrix(mat_C)
  buffer = np.reshape(np.array(buffer), (N, N))
  buffer = buffer + np.transpose(buffer) - np.diag(np.diag(buffer))
  return buffer


### epsilon function
def epsilon(x):
  N = x.shape[0]
  mean = x.mean()
  colmean = np.tile(x.mean(axis = 0), (N, 1))
  rowmean = colmean.T
  return x + mean - colmean - rowmean


### main function
def varGRM_C(trees, file = None, A = math.inf, B = 0, map_func = (lambda x:x), g = (lambda x: 1/(x*(1-x))), var = True):
  if map_func == None: map_func = (lambda x:x)
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
    l = - map_func(tree.interval[0]) + map_func(tree.interval[1])
    if tree.total_branch_length == 0: continue
    for c in tree.nodes():
      descendants = list(tree.samples(c))
      n = len(descendants)
      if(n == 0 or n == N): continue
      t = max(0, min(A, tree.time(tree.parent(c))) - max(B, tree.time(c)))
      mu = l * t * 1e-8
      p = float(n/N)
      matrix.add_square(egrm_C, descendants, mu * g(p))
      if var: matrix.add_square(egrm2_C, descendants, mu * np.power(g(p), 2) * np.power((1 - 2*p), 2))
      if var: tmp1[descendants] += mu * np.power(g(p), 2) * (np.power(p, 2) - 2 * np.power(p, 3))
      if var: tmp2 += mu * np.power(g(p), 2) * np.power(p, 4)
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


def mTMRCA_C(trees, file = None, map_func = (lambda x:x)):
  if map_func == None: map_func = (lambda x:x)
  N = trees.num_samples
  mtmrca_C = matrix.new_matrix(N)
  tmp = 0
  total_l = 0
  pbar = tqdm.tqdm(total = trees.num_trees, 
                   bar_format = '{l_bar}{bar:30}{r_bar}{bar:-30b}',
                   miniters = trees.num_trees // 100,
                   file = file)
  
  for tree in trees.trees():
    l = - map_func(tree.interval[0]) + map_func(tree.interval[1])
    if tree.total_branch_length == 0: continue
    height = 0
    for c in tree.nodes():
      descendants = list(tree.samples(c))
      n = len(descendants)
      if(n == 0 or n == N): continue
      t = tree.time(tree.parent(c)) - tree.time(c)
      height = max(height, tree.time(tree.parent(c)))
      matrix.add_square(mtmrca_C, descendants, t * l)
    tmp += height * l
    total_l += l
    pbar.update(1)
  
  mtmrca = mat_C_to_array(mtmrca_C, N)
  mtmrca = tmp - mtmrca
  mtmrca /= total_l
  pbar.close()
  return mtmrca, total_l


### without C extension (to be added)
def TMRCA(tree):
  N = tree.num_samples()
  rel_nodes = list(tree.nodes())
  times = [tree.time(node) for node in rel_nodes]
  rel_nodes = np.array(rel_nodes)[np.flip(np.argsort(times))]
  tmrca = np.zeros([N, N])
  for c in rel_nodes:
    descendants = list(tree.samples(c))
    n = len(descendants)
    if(n == 0):
      continue
    q = tree.time(c)
    tmrca[np.ix_(descendants, descendants)] = q
  return tmrca

def mTMRCA(trees, file = None):
  N = trees.num_samples
  buffer = np.zeros([N, N])
  total_l = 0
  pbar = tqdm.tqdm(total = trees.num_trees, 
                   bar_format = '{l_bar}{bar:30}{r_bar}{bar:-30b}',
                   miniters = trees.num_trees // 100,
                   file = file)
  for tree in trees.trees():
    if tree.total_branch_length == 0:
      continue
    l = (tree.interval[1] - tree.interval[0])
    total_l += l
    K = TMRCA(tree)
    buffer += K * l
    pbar.update(1)
  buffer /= total_l
  pbar.close()
  return buffer, total_l

def gmTMRCA(trees, file = None):
  N = trees.num_samples
  buffer = np.zeros([N, N])
  total_l = 0
  pbar = tqdm.tqdm(total = trees.num_trees, 
                   bar_format = '{l_bar}{bar:30}{r_bar}{bar:-30b}',
                   miniters = trees.num_trees // 100,
                   file = file)
  for tree in trees.trees():
    if tree.total_branch_length == 0:
      continue
    l = (tree.interval[1] - tree.interval[0])
    total_l += l
    K = TMRCA(tree); np.fill_diagonal(K, 1)
    buffer += np.log(K) * l
    pbar.update(1)
  buffer /= total_l
  buffer = np.exp(buffer); np.fill_diagonal(buffer, 0)
  pbar.close()
  return buffer, total_l

def hmTMRCA(trees, file = None):
  N = trees.num_samples
  buffer = np.zeros([N, N])
  total_l = 0
  pbar = tqdm.tqdm(total = trees.num_trees, 
                   bar_format = '{l_bar}{bar:30}{r_bar}{bar:-30b}',
                   miniters = trees.num_trees // 100,
                   file = file)
  for tree in trees.trees():
    if tree.total_branch_length == 0:
      continue
    l = (tree.interval[1] - tree.interval[0])
    total_l += l
    K = TMRCA(tree); np.fill_diagonal(K, 1)
    buffer += np.power(K, -1) * l
    pbar.update(1)
  buffer /= total_l
  buffer = np.power(buffer, -1); np.fill_diagonal(buffer, 0)
  pbar.close()
  return buffer, total_l
