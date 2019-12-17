### importing packages
import numpy as np

### PLINK
def write_simulation_PLINK(simulation, idx, file):
  X = np.transpose(simulation["hapdata"]["genotype_matrix"][idx])
  y = simulation["phenotypes"]["y"]
  N = X.shape[0]
  
  ped_file = open(file + ".ped", 'a')
  for i in range(N):
    string = "msprime " + str(i+1) + " 0 0 0 " + str(y[i])
    string = string + " ".join(map(str, X[i])) + "\n"
    bytes = ped_file.write(string)
  ped_file.close()
