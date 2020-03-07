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

def test_fun2(N):
  hapdata = simulate_hapdata(l = l, N = N, mutation_rate = mutation_rate, recomb_rate = recomb_rate)
  Z_all = np.transpose(hapdata["genotype_matrix"]).astype("float")
  Z_all -= Z_all.mean(axis=0)
  print("Shape of Z_all: " + str(Z_all.shape))
  print("Number of elements in Z_all: " + str(Z_all.shape[0] * Z_all.shape[1] / 1000000) + "M")
  print("Size of Z_all: " + str(sys.getsizeof(Z_all) / 1000000) + "Mb")
  Z_all /= Z_all.std(axis=0)

for N in (5000, 8000):
  mu = max(memory_usage((test_fun2, (N, ))))
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
  
  for i in idx_trees:
    tree = trees.at_index(i)
    interval = tree.interval
    tl = (interval[1] - interval[0]) * tree.total_branch_length * 1e-8
    total_tl += tl
    EK += getEK(tree) * tl
    pbar.update(1)
  EK /= total_tl
  times.append(datetime.datetime.now())
  return (EK, round(total_tl), times)

def test_fun(N):
  simulation = simulate(N)
  obss = simulation["observations"]["obss"]
  variants = simulation["hapdata"]["variants"]
  flags_obs = get_flags(trees, variants[obss])
  Km, Km_tl, times = getEK_trees(trees, flags_obs)




