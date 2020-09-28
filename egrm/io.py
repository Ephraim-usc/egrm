### importing packages
import os
import math
import numpy as np
import pandas as pd
import tskit
import struct
import tsinfer
import pickle

def getX(hapdata, idx):
  N = hapdata["trees"].num_samples
  X = np.zeros((N, len(idx))).astype("int")
  variants = hapdata["trees"].variants()
  index = 0; num = 0
  for v in variants:
    if index in idx:
      X[:, num] = v.genotypes; num += 1
    index += 1
  return X

### PLINK and GCTA
def write_plink(simulation, samples = None, obss, file):
  from pyplink import PyPlink
  
  N = int(simulation['parameters']['N'] / 2)
  maternals = simulation["diploid"]['maternals']
  paternals = simulation["diploid"]['paternals']
  
  if samples == None:
    samples = np.array(list(range(N)))
  idx_samples = np.concatenate([samples * 2, samples * 2 + 1])
  idx_samples = np.sort(idx_samples)
  
  X = getX(simulation["hapdata"], obss)
  X = X[maternal] + X[paternal]
  X = X[:, idx_samples]
  
  y = simulation["diploid"]["y_diploid"]
  variants = np.ceil(simulation["hapdata"]["loci"][obss])
  
  with PyPlink(file, "w") as bedfile:
    for v in np.transpose(X):
      bedfile.write_genotypes(v)
  
  ped_file = open(file + ".ped", 'a')
  for idx, sample in enumerate(samples):
    string = "FID{}".format(sample+1) + " IID{}".format(sample+1) + " 0 0 0 " + str(y[idx]) + " "
    string = string + " ".join(map(str, X[idx])) + "\n"
    bytes = ped_file.write(string)
  ped_file.close()
  
  fam_file = open(file + ".fam", 'a')
  for idx, sample in enumerate(samples):
    string = "FID{}".format(sample+1) + " IID{}".format(sample+1) + " 0 0 0 " + str(y[idx]) + "\n"
    bytes = fam_file.write(string)
  fam_file.close()
  
  phen_file = open(file + ".phen", 'a')
  for idx, sample in enumerate(samples):
    string = "FID{}".format(sample+1) + " IID{}".format(sample+1) + " " + str(y[idx]) + "\n"
    bytes = phen_file.write(string)
  phen_file.close()
  
  map_file = open(file + ".map", 'a')
  for idx, obs in enumerate(obss):
    string = "1 snp" + str(obs+1) + " 0 " + str(variants[idx]) + "\n"
    bytes = map_file.write(string)
  map_file.close()
  
  bim_file = open(file + ".bim", 'a')
  for idx, obs in enumerate(obss):
    string = "1 snp" + str(obs+1) + " 0 " + str(variants[idx]) + " A T\n"
    bytes = bim_file.write(string)
  bim_file.close()
  
def write_relate(simulation, obss, file): #usually we use obss as idx
  N = int(simulation['parameters']['N'])
  X = getX(simulation["hapdata"], obss)
  variants = np.ceil(simulation["hapdata"]["loci"][obss])
  
  haps_file = open(file + ".haps", 'a')
  i = 0
  for idx, obs in enumerate(obss):
    string = "1 snp" + str(obs+1) + " " + str(variants[idx]) + " A" + " T "
    string = string + " ".join(map(str, X[:, idx])) + "\n"
    bytes = haps_file.write(string)
    i += 1
  haps_file.close()
  
  sample_file = open(file + ".sample",'a')
  sample_file.write("ID_1\tID_2\tmissing\n0\t0\t0\n")
  for idx in range(int(N//2)):
    string = "ind_" + str(idx+1) + "\tind_" + str(idx) + "\t0\n"
    bytes = sample_file.write(string)
  sample_file.close()
  
  map_file = open(file + ".map",'a')
  map_file.write("pos COMBINED_rate Genetic_Map\n")
  for idx, obs in enumerate(obss):
    string = str(variants[idx]) + " " + str(1) + " "
    string = string + str(variants[idx]/1000000) + "\n"
    bytes = map_file.write(string)
  map_file.close()

def write_tsinfer(simulation, obss, file):
  X = getX(simulation["hapdata"], obss)
  variants = np.ceil(simulation["hapdata"]["loci"][obss])
  with tsinfer.SampleData(path=file + ".samples", sequence_length=simulation["parameters"]["l"]) as sample_data:
    for idx, obs in enumerate(obss):
      sample_data.add_site(variants[idx], X[:, idx])
  del X

def read_trees(file):
  return tskit.load(file)

def write_grm(grm, M, file, format = "GCTA"):
  N = grm.shape[0]
  
  if format == "pickle":
    pickle.dump([grm, M], open(file + ".p", "wb" ))
    return
  
  with open("{}.grm.bin".format(file), "wb") as grmfile:
    for idx in range(N):
      for jdx in range(idx + 1):
        val = struct.pack('f', grm[idx, jdx])
        grmfile.write(val)
  
  val = struct.pack('f', int(M))
  with open("{}.grm.N.bin".format(file), "wb") as grmfile:
    for idx in range(N):
      for jdx in range(idx + 1):
        grmfile.write(val)
  
  with open("{}.grm.id".format(file), "w") as grmfile:
    for idx in range(N):
      fid = "FID{}".format(idx)
      iid = "IID{}".format(idx)
      grmfile.write("\t".join([fid, iid]) + os.linesep)

'''
def read_grm(file, format = "GCTA"):
  if format == "pickle":
    grm, M = pickle.load(open(file + ".p", "rb" ))
    return grm, M
  
  BinFileName  = file + ".grm.bin"
  NFileName = file + ".grm.N.bin"
  dt = np.dtype('f4')
  grm_vector = np.fromfile(BinFileName, dtype = dt)
  M_vector = np.fromfile(NFileName, dtype = dt)
  
  N = int(np.floor(np.sqrt(2 * grm_vector.shape[0])))
  grm = np.zeros((N, N))
  M = np.zeros((N, N))
  
  p = 0
  for idx in range(N):
    for jdx in range(idx + 1):
      grm[idx, jdx] = grm_vector[p]
      M[idx, jdx] = M_vector[p]
      p += 1
  
  grm_diags = np.diag(grm); grm = grm.T + grm
  np.fill_diagonal(grm, grm_diags)
  M_diags = np.diag(M); M = M.T + M
  np.fill_diagonal(M, M_diags)
  return (grm, M)
'''
  
'''
def ReadGRMBin(prefix, AllN = False):
    """
    read a GCTA binary GRM file storing relatedness values between individuals
    adapted from an R function on the GCTA website
    Davis McCarthy, February 2015
    """
    BinFileName  = prefix + ".grm.bin"
    NFileName = prefix + ".grm.N.bin"
    IDFileName = prefix + ".grm.id"
    dt = np.dtype('f4') # Relatedness is stored as a float of size 4 in the binary file
    entry_format = 'f' # N is stored as a float in the binary file
    entry_size = calcsize(entry_format)
    ## Read IDs
    ids = pd.read_csv(IDFileName, sep = '\t', header = None)
    ids_vec = ids.iloc[:,1]
    n = len(ids.index)
    ids_diag = ['NA' for x in range(n)]
    n_off = int((n * (n + 1) / 2) - n)
    ids_off = ['NA' for x in range(n_off)]
    ## Generate ids for relatedness values by concatenating individual IDs
    ticker = 0
    for i in range(n):
        for j in range(i):
            if i == j:
                ids_diag[i] = str(ids_vec[i])
            else:
                ids_off[ticker] = str(ids_vec[i]) + '_' + str(ids_vec[j])
                ticker += 1
    ## Read relatedness values
    grm = np.fromfile(BinFileName, dtype = dt)
    ## Read number of markers values
    if AllN:
        N = np.fromfile(NFileName, dtype = dt)
    else:
        with open(NFileName, mode='rb') as f:
            record = f.read(entry_size)
            N = unpack(entry_format, record)[0]
            N = int(N)
    i = sum_n_vec(n)
    out = {'diag': grm[i], 'off': np.delete(grm, i), 'id': ids, 'id_off': ids_off, 'id_diag': ids_diag, 'N': N}
    return(out)
'''

def run_cmd(cmd, log = None):
  if log != None:
    os.system(cmd + " > " + log + " 2>&1")
  else:
    os.system(cmd)

'''
def gcta64(file, log = None):
  run_cmd("/home/rcf-40/caoqifan/cc2/plink-1.07-x86_64/plink --file " + 
            file + " --make-bed --out " + file + " --noweb", log)
  
  run_cmd("/home/rcf-40/caoqifan/cc2/gcta_1.93.0beta/gcta64 --bfile " + 
            file + " --make-grm-bin --out " + file, log)
  
  run_cmd("/home/rcf-40/caoqifan/cc2/gcta_1.93.0beta/gcta64  --reml  --grm " + 
            file + " --out " + file + " --pheno " + file + ".phen", log)

def gcta64reml(file, phen, log = None):
  run_cmd("/home/rcf-40/caoqifan/cc2/gcta_1.93.0beta/gcta64  --reml  --grm " + 
            file + " --out " + file + " --pheno " + phen + ".phen", log)
'''

def cmd_relate(file, log):
  run_cmd("/home1/caoqifan/project/relate_v1.1.2_x86_64_static/bin/Relate --mode All -m 1e-8 -N 30000 --memory 10 --haps " + 
            file + ".haps --sample " + file + ".sample --map " + file + ".map --seed 1 -o " + file, log)
  
  run_cmd("/home1/caoqifan/project/relate_v1.1.2_x86_64_static/bin/RelateFileFormats --mode ConvertToTreeSequence -i " + 
            file + " -o " + file, log)

def cmd_tsinfer(file, log):
  run_cmd("tsinfer infer " + file + ".samples -p -t 4", log)
