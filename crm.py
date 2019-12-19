### import packages
import tskit
import tqdm

### import trees
def read_trees(file):
  trees = tskit.load(file)
  return trees

### coalescent relationship matrix
def g(p):
  return 1/(p*(1-p))

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

def zeta(tree, g):
  N = tree.num_samples()
  rel_nodes = relevant_nodes(tree)
  zetas = np.zeros([N, N])
  for c in rel_nodes:
    descendents = list(tree.samples(c))
    p = len(descendents) / N
    if(p == 0 or p == 1):
      continue
    q = (tree.time(tree.parent(c)) - tree.time(c)) * g(p)
    zetas[np.ix_(descendents, descendents)] += q
  return zetas

def epsilon(x):
  N = x.shape[0]
  mean = x.mean()
  colmean = np.tile(x.mean(axis = 0), (N, 1))
  rowmean = colmean.T
  return x + mean - colmean - rowmean

def getEK(tree):
  buffer = epsilon(zeta(tree, g))
  L = tree.total_branch_length
  return buffer/L

def getEK_trees(trees):
  N = trees.num_samples
  EK = np.zeros([N, N])
  total_length = 0
  pbar = tqdm.tqdm(total = trees.num_trees, 
                   bar_format = '{l_bar}{bar:30}{r_bar}{bar:-30b}',
                   miniters = trees.num_trees // 100)
  for tree in trees.trees():
    interval = tree.interval
    length = interval[1] - interval[0]
    total_length += length
    EK += getEK(tree) * length
    pbar.update(1)
  EK /= total_length
  pbar.close()
  return EK

def filter_trees(trees, variants):
  #needs codes here

