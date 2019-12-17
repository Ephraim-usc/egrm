### importing packages


### PLINK
def write_simulation_PLINK(genotype_matrix, file):
  for i in range(M_obs):
  string = "1 SNP" + str(i+1) + " " + str(math.ceil(variants_obs[i])) + " A" + " T "
  string = string + " ".join(map(str, genotype_matrix_obs[i])) + "\n"
  tmp_file = open(out_file + ".relate/" + out_file + ".haps",'a')
  bytes = tmp_file.write(string)
  tmp_file.close()
