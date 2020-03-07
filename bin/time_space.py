from memory_profiler import memory_usage
from egrm import *
import sys

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





