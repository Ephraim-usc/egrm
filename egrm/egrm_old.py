'''
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

### coalescent relationship matrix
def zeta(tree, A = math.inf, B = 0, g = (lambda x: 1/(x*(1-x)))):
  N = tree.num_samples()
  l = map_func(tree.interval[1]) - map_func(tree.interval[0])
  zetas = np.zeros([N, N])
  total_mu = 0
  for c in tree.nodes():
    descendants = list(tree.samples(c))
    n = len(descendants)
    if(n == 0 or n == N):
      continue
    t = max(0, min(A, tree.time(tree.parent(c))) - max(B, tree.time(c)))
    mu = l * t * 1e-8
    zetas[np.ix_(descendants, descendants)] += mu * g(n / N)
    total_mu += mu
  zetas /= tree.total_branch_length
  return zetas, total_mu

def zeta_C(tree, mat_C, A = math.inf, B = 0, map_func = (lambda x:x), g = (lambda x: 1/(x*(1-x)))):
  N = tree.num_samples()
  l = map_func(tree.interval[1]) - map_func(tree.interval[0])
  total_mu = 0
  for c in tree.nodes():
    descendants = list(tree.samples(c))
    n = len(descendants)
    if (n == 0 or n == N):
      continue
    t = max(0, min(A, tree.time(tree.parent(c))) - max(B, tree.time(c)))
    mu = l * t * 1e-8
    if (t == 0):
      continue
    matrix.add_square(mat_C, descendants, mu * g(n / N))
    total_mu += mu
  return total_mu



def eGRM(trees, file = None, A = math.inf, B = 0, map_func = (lambda x:x), g = (lambda x: 1/(x*(1-x)))):
  if map_func == None:
    map_func = (lambda x:x)
  N = trees.num_samples
  buffer = np.zeros([N, N])
  total_mu = 0
  pbar = tqdm.tqdm(total = trees.num_trees, 
                   bar_format = '{l_bar}{bar:30}{r_bar}{bar:-30b}',
                   miniters = trees.num_trees // 100,
                   file = file)
  for tree in trees.trees():
    if tree.total_branch_length == 0: # especially for trimmed tree sequences
      continue
    K, total_mu_ = zeta(tree, A, B, g)
    total_mu += total_mu_
    buffer += K * total_mu_
    pbar.update(1)
  buffer /= total_mu
  buffer = epsilon(buffer)
  pbar.close()
  return buffer, total_mu

def varGRM(trees, file = None, A = math.inf, B = 0, map_func = (lambda x:x)):
  egrm, egrm_tl = eGRM(trees, file = None, A = A, B = B, map_func = map_func)
  vargrm_, tmp = eGRM(trees, file = None, A = A, B = B, map_func = map_func, g = (lambda x: pow(1/(x*(1-x)), 2)))
  e = np.reciprocal(np.random.poisson(lam=egrm_tl, size=10000).astype("float")).mean()
  vargrm = e * (vargrm_ - np.power(egrm, 2))
  return egrm, vargrm, egrm_tl


def mat_C_to_array(mat_C, N):
  buffer = matrix.export_matrix(mat_C)
  buffer = np.reshape(np.array(buffer), (N, N))
  buffer = buffer + np.transpose(buffer) - np.diag(np.diag(buffer))
  return buffer

def eGRM_C(trees, file = None, A = math.inf, B = 0, map_func = (lambda x:x), g = (lambda x: 1/(x*(1-x)))):
  if map_func == None:
    map_func = (lambda x:x)
  N = trees.num_samples
  total_mu = 0
  mat_C = matrix.new_matrix(N)
  pbar = tqdm.tqdm(total = trees.num_trees, 
                   bar_format = '{l_bar}{bar:30}{r_bar}{bar:-30b}',
                   miniters = trees.num_trees // 100,
                   file = file)
  for tree in trees.trees():
    if tree.total_branch_length == 0: # especially for trimmed tree sequences
      continue
    total_mu += zeta_C(tree, mat_C, A, B, map_func, g)
    pbar.update(1)
  
  buffer = mat_C_to_array(mat_C, N)
  buffer /= total_mu
  pbar.close()
  return buffer, total_mu


##############

def epsilon(x):
  N = x.shape[0]
  mean = x.mean()
  colmean = np.tile(x.mean(axis = 0), (N, 1))
  rowmean = colmean.T
  return x + mean - colmean - rowmean

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





#############

def eGRM_merge(files, N):
  buffer = np.zeros((N, N))
  total_tl = 0
  for file in files:
    with open(file, "rb") as f:
      buffer_, total_tl_ = pickle.load(f)
    buffer += buffer_ * total_tl_
    total_tl += total_tl_
  buffer /= total_tl
  return buffer, total_tl


def _eGRM_C_chunk(trees, index, chunk_size, name):
  N = trees.num_samples
  start = index * chunk_size
  end = min(trees.num_trees, index * chunk_size + chunk_size)
  with open(name + ".log", "a") as f:
    f.write("Computing trees from " + str(start) + " to " + str(end) + "\n")
  
  total_tl = 0
  tree = trees.first()
  mat_C = matrix.new_matrix(N)
  pbar = tqdm.tqdm(total = end - start, 
                   bar_format = '{l_bar}{bar:30}{r_bar}{bar:-30b}',
                   file = open(name + ".log", "a"))
  
  while tree.index < start:
    tree.next()
  while tree.index >= 0 and tree.index < end:
    if tree.total_branch_length == 0: # especially for trimmed tree sequences
      tree.next(); continue
    tl = (tree.interval[1] - tree.interval[0]) * tree.total_branch_length * 1e-8
    total_tl += tl
    zeta_C(tree, mat_C)
    pbar.update(1)
    tree.next()
  
  buffer = mat_C_to_array(mat_C, N)
  matrix.destroy_matrix(mat_C)
  buffer /= total_tl
  buffer = epsilon(buffer)
  pbar.close()
  with open(name + ".log", "a") as f:
    f.write("Done")
  pickle.dump([buffer, total_tl], file = open(name + ".p", "wb"))
  

def eGRM_C_pll(trees, name, cpus = 5):
  N = trees.num_samples
  chunk_size = int(trees.num_trees / cpus) + 1
  print("totally " + str(trees.num_trees) + " trees.")
  os.mkdir(name + ".egrm_tmp"); os.chdir(name + ".egrm_tmp")
  
  ps = []
  for index in range(cpus):
    p = mp.Process(target = _eGRM_C_chunk, args = (trees, index, chunk_size, name + "_" + str(index)))
    ps.append(p); p.start()
    print("New process " + str(p.pid) + " with log file " + name + "_" + str(index) + ".log")
  
  for p in ps:
    p.join()
  
  buffer, total_tl = eGRM_merge([name + "_" + str(index) + ".p" for index in range(cpus)], N)
  
  for file in [name + "_" + str(index) + ".log" for index in range(cpus)] + [name + "_" + str(index) + ".p" for index in range(cpus)]:
    os.remove(file)
  
  os.chdir(".."); os.rmdir(name + ".egrm_tmp")
  print("Done")
  return buffer, total_tl





#####
def eGRM_obs(trees, loci, num = -1):
  N = trees.num_samples
  total_tl = 0
  buffer = np.zeros([N, N])
  tree = trees.first()
  i = 0
  while i < len(loci):
    while tree.interval[0] >= loci[i]:
      i += 1
      if i >= len(loci):
        return buffer/total_tl, total_tl
    while tree.interval[1] < loci[i]:
      tree.next()
    if tree.total_branch_length == 0:
      continue
    tl = (tree.interval[1] - tree.interval[0]) * tree.total_branch_length * 1e-8
    total_tl += tl
    K = zeta(tree)
    buffer += K * tl; #print(str(tree.index) + " " + str(K[0, 0]))
    tree.next(); i += 1
    num -= 1
    if num == 0:
      break
  buffer = epsilon(buffer)
  return buffer/total_tl, total_tl

def TMRCA(tree):
  N = tree.num_samples()
  rel_nodes = list(tree.nodes())
  times = [tree.time(node) for node in rel_nodes]
  rel_nodes = np.array(rel_nodes)[np.flip(np.argsort(times))]
  tmrca = np.zeros([N, N])
  for c in rel_nodes:
    descendants = list(tree.samples(c))
    n = len(descendants)
    if(n == 0 or n == N):
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

def TMRCA_C(tree, mat_C_):
  N = tree.num_samples()
  l = tree.interval[1] - tree.interval[0]
  
  nodes = list(tree.nodes())
  times = [tree.time(node) for node in nodes]
  nodes_sorted = np.array(nodes)[np.flip(np.argsort(times))]
  
  matrix.set_zeros(mat_C_)
  for c in nodes_sorted:
    descendants = list(tree.samples(c))
    n = len(descendants)
    if(n == 0 or n == N):
      continue
    q = tree.time(c)
    matrix.set_square(mat_C_, descendants, float(q * l))
  return l

def mTMRCA_C(trees, file = None):
  N = trees.num_samples
  total_l = 0
  mat_C = matrix.new_matrix(N)
  mat_C_ = matrix.new_matrix(N)
  
  pbar = tqdm.tqdm(total = trees.num_trees, 
                   bar_format = '{l_bar}{bar:30}{r_bar}{bar:-30b}',
                   miniters = trees.num_trees // 100,
                   file = file)
  for tree in trees.trees():
    if tree.total_branch_length == 0:
      continue
    total_l += TMRCA_C(tree, mat_C_)
    matrix.add(mat_C, mat_C_)
    pbar.update(1)
  
  buffer = mat_C_to_array(mat_C, N)
  buffer /= total_l
  pbar.close()
  return buffer, total_l

'''


#########

'''
def relevant_nodes(tree):
  N = tree.num_samples()
  rel_nodes = []
  for i in list(range(N)):
    c = i
    while not c in rel_nodes:
      if c == -1: 
        break
      rel_nodes.append(c)
      c = tree.parent(c)
  rel_nodes.sort()
  return rel_nodes
'''


'''
def affected_nodes(prev_tree, tree):
  buffer = []
  for node in set(prev_tree.nodes()).union(set(tree.nodes())):
    if prev_tree.parent(node) != tree.parent(node):
      buffer.append(node)
    elif sorted(list(prev_tree.samples(node))) != sorted(list(tree.samples(node))):
      buffer.append(node)
  buffer = list(set(buffer))
  return buffer

def zeta_rec(prev_K, tree, prev_tree = None):
  N = tree.num_samples()
  if prev_tree == None:
    tree.prev(); prev_tree = tree.copy(); tree.next()
  K = prev_K * prev_tree.total_branch_length
  affected = affected_nodes(prev_tree, tree)
  
  for node in affected:
    descendants = list(prev_tree.samples(node))
    n = len(descendants)
    if(n == 0 or n == N):
      continue
    time = prev_tree.time(prev_tree.parent(node)) - prev_tree.time(node)
    q = g(n / N)
    K[np.ix_(descendants, descendants)] -= q * time
  
  for node in affected:
    descendants = list(tree.samples(node))
    n = len(descendants)
    if(n == 0 or n == N):
      continue
    time = tree.time(tree.parent(node)) - tree.time(node)
    q = g(n / N)
    K[np.ix_(descendants, descendants)] += q * time
  
  K /= tree.total_branch_length
  return K

def EGRM(trees, num = -1):
  N = trees.num_samples
  buffer = np.zeros([N, N])
  total_tl = 0
  for tree in trees.trees():
    #print(total_tl)
    tl = (tree.interval[1] - tree.interval[0]) * tree.total_branch_length * 1e-8
    total_tl += tl
    K = zeta(tree)
    #print(K[0, 0])
    buffer += K * tl
    num -= 1
    if num == 0:
      break
  buffer /= total_tl
  return buffer, total_tl

def EGRM_rec(trees, num = -1):
  N = trees.num_samples
  buffer = np.zeros([N, N])
  tree = trees.first()
  K = zeta(tree)
  total_tl = (tree.interval[1] - tree.interval[0]) * tree.total_branch_length * 1e-8
  while tree.next():
    tl = (tree.interval[1] - tree.interval[0]) * tree.total_branch_length * 1e-8
    total_tl += tl
    K = zeta_rec(K, tree)
    buffer += K * tl
    num -= 1
    if num == 0:
      break
  buffer /= total_tl
  return buffer, total_tl
'''

'''
def EGRM_obs(trees, loci, num = -1):
  N = trees.num_samples
  total_tl = 0
  buffer = np.zeros([N, N])
  tree = trees.first()
  i = 0
  while i < len(loci):
    while tree.interval[0] >= loci[i]:
      i += 1
      if i >= len(loci):
        return buffer/total_tl, total_tl
    while tree.interval[1] < loci[i]:
      tree.next()
    tl = (tree.interval[1] - tree.interval[0]) * tree.total_branch_length * 1e-8
    total_tl += tl
    K = zeta(tree)
    buffer += K * tl; #print(str(tree.index) + " " + str(K[0, 0]))
    tree.next(); i += 1
    num -= 1
    if num == 0:
      break
  return buffer/total_tl, total_tl

def EGRM_obs_rec(trees, loci, num = -1):
  N = trees.num_samples
  buffer = np.zeros([N, N])
  tree = trees.first()
  i = 0
  while tree.interval[1] < loci[i]:
    tree.next()
  tl = (tree.interval[1] - tree.interval[0]) * tree.total_branch_length * 1e-8
  total_tl = tl
  K = zeta(tree); buffer += K * tl; #print(str(tree.index) + " " + str(K[0, 0]))
  prev_tree = tree.copy(); i += 1
  while i < len(loci):
    while tree.interval[0] >= loci[i]:
      i += 1
      if i >= len(loci):
        return buffer/total_tl, total_tl
    while tree.interval[1] < loci[i]:
      tree.next()
    tl = (tree.interval[1] - tree.interval[0]) * tree.total_branch_length * 1e-8
    total_tl += tl
    K = zeta_rec(K, tree, prev_tree = prev_tree); buffer += K * tl; #print(str(tree.index) + " " + str(K[0, 0]))
    prev_tree = tree.copy()
    tree.next(); i += 1
    num -= 1
    if num == 0:
      break
  return buffer/total_tl, total_tl






def compare(trees1, trees2, n = 100):
  l = min(trees1.last().interval[1], trees2.last().interval[1])
  loci = np.random.randint(0, l, n)
  loci.sort()
  diffs = []
  for locus in loci:
    tree1 = trees1.at(locus)
    tree2 = trees2.at(locus)
    tmrca1 = TMRCA(tree1)
    tmrca2 = TMRCA(tree2)
    tmrca1 /= tmrca1.mean()
    tmrca2 /= tmrca2.mean()
    diff = np.abs(tmrca2 - tmrca1).mean()
    #print("at " + str(locus) + " diff " + str(diff))
    diffs.append(diff)
  return diffs 
'''
