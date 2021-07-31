### import packages
import tskit
import tqdm
import math
import numpy as np
import pandas as pd


### import C extension
import matrix


# exports [matrix._matrix_C_API] object into a 2d numpy array
# mat_C: [matrix._matrix_C_API] object initiated by matrix.new_matrix, only stores the upper triangle elements of a square matrix
# N: the number of columns/rows
def mat_C_to_ndarray(mat_C, N): 
  buffer = matrix.export_ndarray(mat_C)
  buffer = buffer + np.transpose(buffer) - np.diag(np.diag(buffer))
  return buffer


### defines [Gmap] object which maps the physical position (in bp) into genetic position (in unit of 10^-6 cM)
### can be initiated by Gmap(filename), where filename is a (comma/space/tab separated) three-column file 
### with first column specifying the physical position in bp and the third column specifying the genetic position in cM. 
### The second column is not used. The first line will always be ignored as the header.
class Gmap:
  def __init__(self, filename):
    if filename is None:
      self.mapped = False
      return
    self.table = pd.read_csv(filename, sep = None, engine = 'python')
    self.pos = self.table.iloc[:, 0].astype(int)
    self.gpos = self.table.iloc[:, 2].astype(float) * 1e6
    self.max = self.table.shape[0]
    self.i = 0
    self.mapped = True
  
  def __call__(self, pos):
    if self.mapped == False:
      return pos
    while (self.i > 0 and pos < self.pos[self.i - 1]):
      self.i -= 1
    while (self.i < self.max and pos > self.pos[self.i]):
      self.i += 1
    if self.i == 0:
      return 0
    if self.i >= self.max:
      return self.gpos[self.max - 1]
    A = pos - self.pos[self.i-1]
    B = (self.gpos[self.i] - self.gpos[self.i-1])/(self.pos[self.i] - self.pos[self.i-1])
    C = self.gpos[self.i-1]
    return A * B + C


### main functions

# computes the eGRM (and varGRM if var == True) based on the tskit TreeSequence trees.
# trees: tskit TreeSequence object.
# log: tqdm log file path
# rlim, alim: most recent and most ancient time limits (in unit of generation) between which the eGRM (varGRM) is computed.
# left, right: leftmost and rightmost positions (in unit of base pair) between which the eGRM (varGRM) is computed.
# gmap: [Gmap] object that maps the physical position (in unit of base pair) into genetic position (in unit of 10^-6 centimorgan).
# g: a scaling function of the allele frequency. Use the default for the standard GRM and eGRM definitions.
# var: True/False variable indicating whether the varGRM is to be computed.
# sft: True/False variable indicating whether the first tree is to be skipped. Not recommended to use together with [left] and [right] options.
def varGRM_C(trees, log = None, 
             rlim = 0, alim = math.inf, 
             left = 0, right = math.inf, 
             gmap = Gmap(None), g = (lambda x: 1/(x*(1-x))), 
             var = True, sft = False):
  
  N = trees.num_samples
  egrm_C = matrix.new_matrix(N)
  if var:
    egrm2_C = matrix.new_matrix(N)
    tmp1 = np.zeros(N)
    tmp2 = 0
  
  total_mu = 0
  pbar = tqdm.tqdm(total = trees.num_trees, 
                   bar_format = '{l_bar}{bar:30}{r_bar}{bar:-30b}',
                   miniters = trees.num_trees // 100,
                   file = log)
  
  trees_ = trees.trees()
  if sft: next(trees_)
  
  for tree in trees_:
    pbar.update(1)
    if tree.total_branch_length == 0: continue
    l = - gmap(max(left, tree.interval[0])) + gmap(min(right, tree.interval[1]))
    if l <= 0: continue
    
    for c in tree.nodes():
      descendants = list(tree.samples(c))
      n = len(descendants)
      if(n == 0 or n == N): continue
      t = max(0, min(alim, tree.time(tree.parent(c))) - max(rlim, tree.time(c)))
      mu = l * t * 1e-8
      p = float(n/N)
      matrix.add_square(egrm_C, descendants, mu * g(p))
      if var:
        matrix.add_square(egrm2_C, descendants, mu * np.power(g(p), 2) * np.power((1 - 2*p), 2))
        tmp1[descendants] += mu * np.power(g(p), 2) * (np.power(p, 2) - 2 * np.power(p, 3))
        tmp2 += mu * np.power(g(p), 2) * np.power(p, 4)
      total_mu += mu
  
  egrm = mat_C_to_ndarray(egrm_C, N)
  if var:
    egrm2 = mat_C_to_ndarray(egrm2_C, N) 
  
  egrm /= total_mu
  if var:
    egrm2 /= total_mu
    tmp1 /= total_mu
    tmp2 /= total_mu
  
  egrm -= egrm.mean(axis = 0)
  egrm -= egrm.mean(axis = 1, keepdims=True)
  if var:
    e = np.reciprocal((lambda x:x[x!=0])(np.random.poisson(lam=total_mu, size=10000)).astype("float")).mean()
    vargrm = e * (egrm2 + np.tile(tmp1, (N, 1)) + np.tile(tmp1, (N, 1)).transpose() + tmp2 - np.power(egrm, 2))
  else:
    vargrm = None
  
  pbar.close()
  return egrm, vargrm, total_mu


# computes the mean TMRCA (mTMRCA) based on the tskit TreeSequence [trees].
# trees: tskit TreeSequence object.
# log: tqdm log file path
# left, right: leftmost and rightmost positions (in unit of base pair) between which the mTMRCA is computed.
# gmap: Gmap object that maps the physical position (in unit of base pair) into genetic position (in unit of 10^-6 centimorgan).
# sft: True/False variable indicating whether the first tree is to be skipped. Not recommended to use together with [left] and [right] options.
def mTMRCA_C(trees, log = None, 
             left = 0, right = math.inf, 
             gmap = Gmap(None), sft = False):
  
  N = trees.num_samples
  mtmrca_C = matrix.new_matrix(N)
  tmp = 0
  total_l = 0
  pbar = tqdm.tqdm(total = trees.num_trees, 
                   bar_format = '{l_bar}{bar:30}{r_bar}{bar:-30b}',
                   miniters = trees.num_trees // 100,
                   file = log)
  
  trees_ = trees.trees()
  if sft: next(trees_)
  
  for tree in trees_:
    pbar.update(1)
    if tree.total_branch_length == 0: continue
    l = - gmap(max(left, tree.interval[0])) + gmap(min(right, tree.interval[1]))
    if l <= 0: continue
    
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
  
  mtmrca = mat_C_to_ndarray(mtmrca_C, N)
  mtmrca = tmp - mtmrca
  mtmrca /= total_l
  pbar.close()
  return mtmrca, total_l


### without C extension

# the non-C version of varGRM_C
def varGRM(trees, log = None, 
             rlim = 0, alim = math.inf, 
             left = 0, right = math.inf, 
             gmap = Gmap(None), g = (lambda x: 1/(x*(1-x))), 
             var = True, sft = False):
  
  N = trees.num_samples
  egrm = np.zeros([N, N])
  if var:
    egrm2 = np.zeros([N, N])
    tmp1 = np.zeros(N)
    tmp2 = 0
  
  total_mu = 0
  pbar = tqdm.tqdm(total = trees.num_trees, 
                   bar_format = '{l_bar}{bar:30}{r_bar}{bar:-30b}',
                   miniters = trees.num_trees // 100,
                   file = log)
  
  trees_ = trees.trees()
  if sft: next(trees_)
  
  for tree in trees_:
    pbar.update(1)
    if tree.total_branch_length == 0: continue
    l = - gmap(max(left, tree.interval[0])) + gmap(min(right, tree.interval[1]))
    if l <= 0: continue
    
    for c in tree.nodes():
      descendants = list(tree.samples(c))
      n = len(descendants)
      if(n == 0 or n == N): continue
      t = max(0, min(alim, tree.time(tree.parent(c))) - max(rlim, tree.time(c)))
      mu = l * t * 1e-8
      p = float(n/N)
      egrm[np.ix_(descendants, descendants)] += mu * g(p)
      if var:
        egrm2[np.ix_(descendants, descendants)] += mu * np.power(g(p), 2) * np.power((1 - 2*p), 2)
        tmp1[descendants] += mu * np.power(g(p), 2) * (np.power(p, 2) - 2 * np.power(p, 3))
        tmp2 += mu * np.power(g(p), 2) * np.power(p, 4)
      total_mu += mu
  
  egrm /= total_mu
  if var:
    egrm2 /= total_mu
    tmp1 /= total_mu
    tmp2 /= total_mu
  
  egrm -= egrm.mean(axis = 0)
  egrm -= egrm.mean(axis = 1, keepdims=True)
  if var:
    e = np.reciprocal((lambda x:x[x!=0])(np.random.poisson(lam=total_mu, size=10000)).astype("float")).mean()
    vargrm = e * (egrm2 + np.tile(tmp1, (N, 1)) + np.tile(tmp1, (N, 1)).transpose() + tmp2 - np.power(egrm, 2))
  else:
    vargrm = None
  
  pbar.close()
  return egrm, vargrm, total_mu


# the non-C version of mTMRCA_C
def mTMRCA(trees, log = None, 
             left = 0, right = math.inf, 
             gmap = Gmap(None), sft = False):
  
  N = trees.num_samples
  mtmrca = np.zeros([N, N])
  tmp = 0
  total_l = 0
  pbar = tqdm.tqdm(total = trees.num_trees, 
                   bar_format = '{l_bar}{bar:30}{r_bar}{bar:-30b}',
                   miniters = trees.num_trees // 100,
                   file = log)
  
  trees_ = trees.trees()
  if sft: next(trees_)
  
  for tree in trees_:
    pbar.update(1)
    if tree.total_branch_length == 0: continue
    l = - gmap(max(left, tree.interval[0])) + gmap(min(right, tree.interval[1]))
    if l <= 0: continue
    
    height = 0
    for c in tree.nodes():
      descendants = list(tree.samples(c))
      n = len(descendants)
      if(n == 0 or n == N): continue
      t = tree.time(tree.parent(c)) - tree.time(c)
      height = max(height, tree.time(tree.parent(c)))
      mtmrca[np.ix_(descendants, descendants)] += t * l
    tmp += height * l
    total_l += l
  
  mtmrca = tmp - mtmrca
  mtmrca /= total_l
  pbar.close()
  return mtmrca, total_l
