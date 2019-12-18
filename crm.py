### import packages
import tskit

### import trees
def read_trees(file):
  trees = tskit.load(file)
  return trees

### coalescent relationship matrix
def g(p):
  return 1/(p*(1-p))

def relevant_nodes(tree):
  N = tree.num_samples
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
  N = tree.num_samples
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
  mean = x.mean()
  colmean = np.tile(x.mean(axis = 0), (N, 1))
  rowmean = tmp2.T
  return x + mean - colmean - rowmean

def getEK(tree):
  buffer = epsilon(zeta(tree, g))
  L = tree.total_branch_length
  return buffer/L

def getEK_trees(trees): #scale的问题依然需要修改
  EK = np.zeros([N, N])
  cnt = 0
  for tree in trees.trees():
    cnt += 1
    if cnt % 100 == 0:
      printf("tree " + str(cnt) + "/" + str(trees.num_trees))
    
    interval = tree.interval
    obs = False
    for v in variants_obs:
      if interval[0] < v < interval[1]:
        obs = True
    if obs == False:
      continue
    
    length = interval[1] - interval[0]
    EK += getEK(tree) * length
  return(EK)
