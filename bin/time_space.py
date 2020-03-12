from memory_profiler import memory_usage
from egrm import *
import sys
import datetime

def test_fun1(N):
  hapdata = simulate_hapdata(l = l, N = N, mutation_rate = mutation_rate, recomb_rate = recomb_rate)
  Z_all = np.transpose(hapdata["genotype_matrix"]).astype("float")
  Z_all -= Z_all.mean(axis=0)
  print("Shape of Z_all: " + str(Z_all.shape))
  print("Number of elements in Z_all: " + str(Z_all.shape[0] * Z_all.shape[1] / 1000000) + "M")
  print("Size of Z_all: " + str(sys.getsizeof(Z_all) / 1000000) + "Mb")
  #Z_all /= Z_all.std(axis=0)

def test_fun(N):
  hapdata = simulate_hapdata(l = l, N = N, mutation_rate = mutation_rate, recomb_rate = recomb_rate)
  print("N: " + str(N))
  print("M: " + str(hapdata["M"]))
  print("N * M: " + str(hapdata["M"] * N / 1000000) + "M")
  Z_all = np.transpose(hapdata["genotype_matrix"]).astype("float")
  Z_all -= Z_all.mean(axis=0)
  Z_all /= Z_all.std(axis=0)
  print("Shape of Z_all: " + str(Z_all.shape))
  K_all = np.dot(Z_all, np.transpose(Z_all)) / Z_all.shape[1]

for N in (1000, 2000, 5000, 10000):
  mu = max(memory_usage((test_fun, (N, ))))
  print("Memory usage: " + str(mu) + "Mb")
  print(" ")


def getEK_trees(trees, flags = None, file = None):
  if (flags == None):
    flags = [True] * trees.num_trees
  elif len(flags) != trees.num_trees:
    print("incorrect flags length!")
  idx_trees = np.where(flags)[0].tolist()
  
  N = trees.num_samples
  EK = np.zeros([N, N])
  total_tl = 0
  times = []
  
  times.append(datetime.datetime.now())
  for i in idx_trees:
    tree = trees.at_index(i)
    interval = tree.interval
    tl = (interval[1] - interval[0]) * tree.total_branch_length * 1e-8
    total_tl += tl
    EK += getEK(tree) * tl
    info = str(i) + "/" + str(trees.num_trees) + "   " + str(datetime.datetime.now() - times[-1])
    print(info)
    times.append(datetime.datetime.now())
  EK /= total_tl
  return (EK, round(total_tl), times)



def test_fun(N):
  simulation = simulate(N = N)
  trees = simulation["hapdata"]["tree_sequence"]
  #obss = simulation["observations"]["obss"]
  obss = np.array(list(range(100)))
  variants = simulation["hapdata"]["variants"]
  flags_obs = get_flags(trees, variants[obss])
  Km, Km_tl, times = getEK_trees(trees, flags_obs)
  return times

N = 10000
simulation = simulate(N = N)
trees = simulation["hapdata"]["tree_sequence"]

def test_fun(n):
  obss = np.array(list(range(n)))
  print(obss)
  variants = simulation["hapdata"]["variants"]
  flags_obs = get_flags(trees, variants[obss])
  print("N: " + str(N))
  Km, Km_tl, times = getEK_trees(trees, flags_obs)

for n in (100,):
  print("Number of trees to process: " + str(n))
  mu = max(memory_usage((test_fun, (n, ))))
  print("Memory usage: " + str(mu) + "Mb")
  print(" ")







lengths = []
for tree in trees.trees():
  length = tree.interval[1] - tree.interval[0]
  lengths.append(length)

lengths.sort(reverse=True)
lengths = np.array(lengths)



