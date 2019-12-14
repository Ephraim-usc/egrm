### importing packages
import sys
import time
import datetime
import msprime
import tskit
import itertools
import math
import numpy as np
import pandas as pd
import os

### fixed parameters
l = 32000000
N = 1000
mutation_rate = 1e-8
recomb_rate = 1e-8

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
if sys.argv == ['']:
  out_file = "0.1_0.2_-0.5_1_1_rep_tmp"
  cas_ratio = 0.1
  obs_ratio = 0.2
  Alpha = -0.5
  Beta = 1
  h2g = 1.0
  def printf(string):
    print(string)
else:
  out_file = str(sys.argv[1])
  cas_ratio = float(sys.argv[2])
  obs_ratio = float(sys.argv[3])
  Alpha = float(sys.argv[4])
  Beta = float(sys.argv[5])
  h2g = float(sys.argv[6])
  def printf(string):
    out = open(out_file, "a")
    out.write(str(string) + "\n")
    out.close()

### msprime
population_configurations = [msprime.PopulationConfiguration(sample_size = N, initial_size = N_EU, growth_rate = r_EU)]

demo_events = [msprime.PopulationParametersChange(time=T_EU_AS, initial_size = N_B, growth_rate=0),
              msprime.PopulationParametersChange(time=T_B, initial_size = N_AF, growth_rate=0),
              msprime.PopulationParametersChange(time=T_AF, initial_size = N_A, growth_rate=0)]

tree_sequence = msprime.simulate(length = l, population_configurations = population_configurations,
                recombination_rate = recomb_rate, mutation_rate = mutation_rate, demographic_events = demo_events)

genotype_matrix = tree_sequence.genotype_matrix().astype("uint16")
MAFs = genotype_matrix.mean(axis = 1); MACs = np.array(list(map(min, MAFs, 1-MAFs)))
variants = np.array([v.position for v in tree_sequence.variants()])
M = len(variants)

### phenotype simulation
M_cas = int(round(M * cas_ratio))
weights = np.multiply(np.power(MAFs, Alpha), np.power(1-MAFs, Alpha))
weights /= weights.sum()
cass = np.random.choice(range(M), M_cas, replace = False, p = weights); cass.sort()

X_cas = np.transpose(genotype_matrix[cass])
Z_cas = X_cas - np.mean(X_cas, axis=0)
Z_cas = Z_cas / np.std(Z_cas, axis=0)

betas = np.random.normal(scale = np.sqrt(h2g / M_cas), size = (M_cas, 1))

G = np.dot(X_cas, betas)
G.shape = (N, )
s2g = np.var(G, ddof = 1)
s2e = s2g * (1/h2g - 1)

e = np.random.normal(scale=np.sqrt(s2e), size=N)
y = G + e

### observed data
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

genotype_matrix_obs = genotype_matrix[obss]
variants_obs = variants[obss]

### output maf
MAFs_obs = MAFs[obss]
hist, bin_edges = np.histogram(MAFs, np.arange(0, 0.51, 0.01), density = False)
hist_5, bin_edges = np.histogram(MAFs_5, np.arange(0, 0.51, 0.01), density = False)
hist_obs, bin_edges = np.histogram(MAFs_obs, np.arange(0, 0.51, 0.01), density = False)

MAFs_cas = MAFs[cass]
hist_cas, bin_edges = np.histogram(MAFs_cas, np.arange(0, 0.51, 0.01), density = False)

hists = pd.DataFrame({"total":hist, "ge5":hist_5, "obs":hist_obs, "cas":hist_cas})
f = open(out_file + ".maf", 'a')
f.write(hists.to_string())
f.close()

### relate
os.system("mkdir " + out_file + ".relate")
printf("Writing relate input files ...")

for i in range(M_obs):
  string = "1 SNP" + str(i+1) + " " + str(math.ceil(variants_obs[i])) + " A" + " T "
  string = string + " ".join(map(str, genotype_matrix_obs[i])) + "\n"
  tmp_file = open(out_file + ".relate/" + out_file + ".haps",'a')
  bytes = tmp_file.write(string)
  tmp_file.close()

tmp_file = open(out_file + ".relate/" + out_file + ".sample",'a')
tmp_file.write("ID_1\tID_2\tmissing\n")
tmp_file.write("0\t0\t0\n")
tmp_file.close()
for i in range(int(N/2)):
  string = "msp_" + str(i) + "\tmsp_" + str(i) + "\t0\n"
  tmp_file = open(out_file + ".relate/" + out_file + ".sample",'a')
  bytes = tmp_file.write(string)
  tmp_file.close()

tmp_file = open(out_file + ".relate/" + out_file + ".map",'a')
tmp_file.write("pos COMBINED_rate Genetic_Map\n")
tmp_file.close()
for i in range(M_obs):
  string = str(math.ceil(variants_obs[i])) + " " + str(1) + " " + \
           str(math.ceil(variants_obs[i])/1000000) + "\n"
  tmp_file = open(out_file + ".relate/" + out_file + ".map",'a')
  bytes = tmp_file.write(string)
  tmp_file.close()

relate_cmd = "cd " + out_file + ".relate; " + "/home/rcf-40/caoqifan/cc2/relate_v1.0.16_x86_64_static/bin/Relate --mode All -m 1e-8 -N 30000 --haps " + \
out_file + ".haps --sample " + out_file + ".sample --map " + out_file + ".map --seed 1 -o " + out_file
os.system(relate_cmd)

relate_format = "/home/rcf-40/caoqifan/cc2/relate_v1.0.16_x86_64_static/bin/RelateFileFormats --mode ConvertToTreeSequence -i " + \
out_file + ".relate/" + out_file + " -o " + out_file + ".relate/" + out_file
os.system(relate_format)

relate_ts = tskit.load(out_file + ".relate/" + out_file + ".trees")

### tsinfer
import tsinfer

sample_data = tsinfer.SampleData(sequence_length=l)

for i in range(M_obs):
  sample_data.add_site(variants_obs[i], X_obs[:, i])

tsinfer_ts = tsinfer.infer(sample_data)

### K matrices
def g(p):
  return 1/(p*(1-p))

def relevant_nodes(tree):
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


a = time.time() - time.time()
b = a
def zeta_(tree, g):
  global a; global b
  a -= time.time()
  
  num_nodes = tree.num_nodes
  times = np.zeros(num_nodes+1)
  times[-1] = tree.time(tree.root)
  num_descendants = np.zeros(num_nodes)
  for c in range(num_nodes):
    times[c] = tree.time(c)
    b -= time.time()
    num_descendants[c] = len(list(tree.samples(c)))
    b += time.time()
  relevant_nodes = np.where(num_descendants > 0)[0]
  
  zetas = np.zeros([N, N])
  ps = np.zeros(num_nodes)
  qs = np.zeros(num_nodes)
  for c in relevant_nodes:
    descendents = list(tree.samples(c))
    ps[c] = len(descendents) / N
    if(ps[c] == 0 or ps[c] == 1):
      continue
    qs[c] = (times[tree.parent(c)] - times[c]) * g(ps[c])
    zetas[np.ix_(descendents, descendents)] += qs[c]
  a += time.time()
  print("a: " + str(a))
  print("b: " + str(b))
  return zetas

def epsilon(x):
  tmp1 = x.mean()
  tmp2 = np.tile(x.mean(axis = 0), (N, 1))
  tmp3 = tmp2.T
  buffer = x + tmp1 - tmp2 - tmp3
  return buffer

def getEK(tree):
  buffer = epsilon(zeta(tree, g))
  L = tree.total_branch_length
  return buffer/L

def getEK_trees(trees):
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

Km = getEK_trees(tree_sequence)
Km_tsinfer = getEK_trees(tsinfer_ts)
Km_relate = getEK_trees(relate_ts)

### prep for prediction
K_cas = np.dot(X_cas, np.transpose(X_cas))
K_obs = np.dot(X_obs, np.transpose(X_obs))

diags = np.diag_indices(N)
non_diags = np.where(~np.eye(N,dtype=bool))

table = {"K_cas":K_cas[non_diags].flatten(), "K_obs":K_obs[non_diags].flatten(), \
"Km":Km[non_diags].flatten(), "Km_tsinfer":Km_tsinfer[non_diags].flatten(), "Km_relate":Km_relate[non_diags].flatten()}

table = pd.DataFrame(data=table)
printf(table.corr(method ='pearson'))

###
def BLUP(K, y_train, trains, tests, M_, h2 = 0.9):
  N_train = len(trains)
  
  I = (1/h2 - h2) * M_ * np.identity(N_train)
  V = K[trains, trains] + I
  
  Kt = K[tests, :][:, trains]
  
  Vi = np.linalg.inv(V)
  y_ = np.dot(Kt, np.dot(Vi, y_train))
  return y_

a = []
b = []
c = []
d = []
e = []
for i in range(20):
  tests = np.random.choice(N, math.floor(N * 0.25), replace = False)
  tests.sort()
  trains = [i for i in range(N) if i not in tests]
  
  y_train = y[trains]
  y_test = y[tests]
  
  y_ = BLUP(K_cas, y_train, trains, tests, M_cas, h2 = 0.9)
  a.append(np.corrcoef(y_, y_test)[0, 1])
  
  y_ = BLUP(K_obs, y_train, trains, tests, M_obs, h2 = 0.9)
  b.append(np.corrcoef(y_, y_test)[0, 1])
  
  y_ = BLUP(Km, y_train, trains, tests, M_obs, h2 = 0.9)
  c.append(np.corrcoef(y_, y_test)[0, 1])
  
  y_ = BLUP(Km_tsinfer, y_train, trains, tests, M_obs, h2 = 0.9)
  d.append(np.corrcoef(y_, y_test)[0, 1])
  
  y_ = BLUP(Km_relate, y_train, trains, tests, M_obs, h2 = 0.9)
  e.append(np.corrcoef(y_, y_test)[0, 1])


###### report
a = np.array(a)
b = np.array(b)
c = np.array(c)
d = np.array(d)
e = np.array(e)

printf(" ")
printf("l = " + str(l))
printf("N = " + str(N))
printf("M = " + str(M))
printf("M_cas / M = " + str(np.array(M_cas/M).round(2))) #输出部分急需修改！！！
printf("M_obs / M_5 = " + str(np.array(M_obs/M_5).round(2)))
printf("h2g = " + str(h2g))

printf("true number of trees = " + str(tree_sequence.num_trees))
printf("tsinfer inferred number of trees = " + str(tsinfer_ts.num_trees))
printf("relate inferred number of trees = " + str(relate_ts.num_trees))

printf("K_cas BLUP: " + str(a.mean().round(3)) + " +- " + str(a.std().round(3)))
printf("K_obs BLUP: " + str(b.mean().round(3)) + " +- " + str(b.std().round(3)))
printf("Km BLUP: " + str(c.mean().round(3)) + " +- " + str(c.std().round(3)))
printf("Km_tsinfer BLUP: " + str(d.mean().round(3)) + " +- " + str(d.std().round(3)))
printf("Km_relate BLUP: " + str(e.mean().round(3)) + " +- " + str(e.std().round(3)))
