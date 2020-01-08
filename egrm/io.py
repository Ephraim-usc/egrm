### importing packages
import os
import math
import numpy as np
import pandas as pd
import tskit
import struct

### PLINK and GCTA
def write_plink(simulation, idx, file):
  idx_diploid = idx * 2; idx_diploid = np.concatenate([idx_diploid, idx_diploid + 1])
  idx_diploid = np.sort(idx_diploid)
  X = np.transpose(simulation["hapdata"]["genotype_matrix_diploid"][idx_diploid])
  y = simulation["phenotypes"]["y_diploid"]
  variants = simulation["hapdata"]["variants"]
  N = X.shape[0]
  
  ped_file = open(file + ".ped", 'a')
  for i in range(N):
    string = "FID{}".format(i+1) + "IID{}".format(i+1) + " 0 0 0 " + str(y[i]) + " "
    string = string + " ".join(map(str, X[i])) + "\n"
    bytes = ped_file.write(string)
  ped_file.close()
  
  fam_file = open(file + ".fam", 'a')
  for i in range(N):
    string = "FID{}".format(i+1) + "IID{}".format(i+1) + " 0 0 0 " + str(y[i]) + "\n"
    bytes = fam_file.write(string)
  fam_file.close()
  
  phen_file = open(file + ".phen", 'a')
  for i in range(N):
    string = "FID{}".format(i+1) + "IID{}".format(i+1) + " " + str(y[i]) + "\n"
    bytes = phen_file.write(string)
  phen_file.close()
  
  map_file = open(file + ".map", 'a')
  for i in idx:
    string = "1 snp" + str(i+1) + " 0 " + str(variants[i]) + "\n"
    bytes = map_file.write(string)
  map_file.close()
  
  bim_file = open(file + ".bim", 'a')
  for i in idx:
    string = "1 snp" + str(i+1) + " 0 " + str(variants[i]) + " A T\n"
    bytes = bim_file.write(string)
  bim_file.close()
  
def write_relate(simulation, idx, file): #usually we use obss as idx
  M = len(idx)
  genotype_matrix = simulation["hapdata"]["genotype_matrix"]
  variants = simulation["hapdata"]["variants"]
  N = genotype_matrix.shape[1]
  
  haps_file = open(file + ".haps", 'a')
  for i in idx:
    string = "1 snp" + str(i+1) + " " + str(math.ceil(variants[i])) + " A" + " T "
    string = string + " ".join(map(str, genotype_matrix[i])) + "\n"
    bytes = haps_file.write(string)
  haps_file.close()
  
  sample_file = open(file + ".sample",'a')
  sample_file.write("ID_1\tID_2\tmissing\n0\t0\t0\n")
  for i in range(int(N//2)):
    string = "ind_" + str(i+1) + "\tind_" + str(i) + "\t0\n"
    bytes = sample_file.write(string)
  sample_file.close()
  
  map_file = open(file + ".map",'a')
  map_file.write("pos COMBINED_rate Genetic_Map\n")
  for i in idx:
    string = str(math.ceil(variants[i])) + " " + str(1) + " "
    string = string + str(math.ceil(variants[i])/1000000) + "\n"
    bytes = map_file.write(string)
  map_file.close()

def read_trees(file):
  return tskit.load(file)

def write_grm(grm, M, file):
  N = grm.shape[0]
  with open("{}.grm.bin".format(file), "wb") as grmfile:
    for idx in range(N):
      for jdx in range(idx + 1):
        val = struct.pack('f', grm[idx, jdx])
        grmfile.write(val)
  
  val = struct.pack('i', int(M))
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
def write_crm(grm, file):
  N = crm.shape[0]
  crm_df = pd.DataFrame(data = crm)
  crm_df["row"] = range(N)
  crm_melted = pd.melt(crm_df, id_vars=['row'], var_name='column', value_name='crm')
  
  crm_file = open(file + ".crm", 'a')
  for i in range(N):
    vector = crm_melted.iloc[i, :]
    string = str(vector[0] + 1) + " " + str(vector[1] + 1) + " 1000000 " + str(vector[2]) + "\n"
    bytes = crm_file.write(string)
  crm_file.close()
'''

def run_cmd(cmd, log = None):
  if log != None:
    os.system(cmd + " > " + log + " 2>&1")
  else:
    os.system(cmd)

def gcta64(file, log = None):
  run_cmd("/home/rcf-40/caoqifan/cc2/plink-1.07-x86_64/plink --file " + 
            file + " --make-bed --out " + file + " --noweb", log)
  
  run_cmd("/home/rcf-40/caoqifan/cc2/gcta_1.93.0beta/gcta64 --bfile " + 
            file + " --make-grm-bin --out " + file, log)
  
  run_cmd("/home/rcf-40/caoqifan/cc2/gcta_1.93.0beta/gcta64  --reml  --grm " + 
            file + " --out " + file + " --pheno " + file + ".phen", log)
  
def relate(file, log):
  run_cmd("/home/rcf-40/caoqifan/cc2/relate_v1.0.16_x86_64_static/bin/Relate --mode All -m 1e-8 -N 30000 --haps " + 
            file + ".haps --sample " + file + ".sample --map " + file + ".map --seed 1 -o " + file, log)
  
  run_cmd("/home/rcf-40/caoqifan/cc2/relate_v1.0.16_x86_64_static/bin/RelateFileFormats --mode ConvertToTreeSequence -i " + 
            file + " -o " + file, log)
