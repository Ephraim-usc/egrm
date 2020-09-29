import numpy as np
import pandas as pd
import math

def BLUP(K, y_train, trains, tests, h2 = 0.9):
  N_train = len(trains)
  
  I = (1/h2 - h2) * np.identity(N_train)
  V = K[trains, :][:, trains] + I
  
  Kt = K[tests, :][:, trains]
  
  Vi = np.linalg.inv(V)
  y_ = np.dot(Kt, np.dot(Vi, y_train))
  return y_

def getKs(simulation):
  Ks = {}
  if "K_all" in simulation["Ks"]:
    Ks["K_all"] = simulation["Ks"]["K_all"]
  if "K_cas" in simulation["Ks"]:
    Ks["K_cas"] = simulation["Ks"]["K_cas"]
  if "K_obs" in simulation["Ks"]:
    Ks["K_obs"] = simulation["Ks"]["K_obs"]
  if "EK" in simulation["Ks"]:
    Ks["EK"] = simulation["Ks"]["EK"]
  if "EK_relate" in simulation["Ks"]:
    Ks["EK_relate"] = simulation["Ks"]["EK_relate"]
  if "EK_relate_phased" in simulation["Ks"]:
    Ks["EK_relate_phased"] = simulation["Ks"]["EK_relate_phased"]
  if "EK_tsinfer" in simulation["Ks"]:
    Ks["EK_tsinfer"] = simulation["Ks"]["EK_tsinfer"]
  if "EK_tsdate" in simulation["Ks"]:
    Ks["EK_tsdate"] = simulation["Ks"]["EK_tsdate"]
  if "mTMRCA" in simulation["Ks"]:
    Ks["mTMRCA"] = simulation["Ks"]["mTMRCA"]
  return Ks

def phenotype_impute(simulation, repeats = 100):
  y = simulation['phenotypes']['y']
  N = simulation["parameters"]["N"]
  Ks = getKs(simulation)
  
  results = {}
  for key in Ks.keys():
      results[key] = []
  
  for i in range(repeats):
    h2g = simulation["parameters"]["h2g"]
    tests = np.random.choice(N, math.floor(N * 0.25), replace = False)
    tests.sort()
    trains = np.array([i for i in range(N) if i not in tests])
    y_train = y[trains]
    y_test = y[tests]
    for key in Ks.keys():
      y_ = BLUP(Ks[key], y_train, trains, tests, h2 = h2g)
      results[key].append(np.corrcoef(y_, y_test)[0, 1])
  
  for key in Ks.keys():
      results[key] = np.array(results[key])
  
  return results

def K_diploid(K, maternals, paternals):
  K1 = K[maternals, :][:, maternals]
  K2 = K[maternals, :][:, paternals]
  K3 = K[paternals, :][:, maternals]
  K4 = K[paternals, :][:, paternals]
  return (K1 + K2 + K3 + K4)/2

'''
def phenotype_impute_diploid(simulation, repeats = 1000):
  N = simulation["parameters"]["N"]
  y_diploid = simulation['diploid']['y_diploid']
  maternals = simulation['diploid']['maternals']
  paternals = simulation['diploid']['paternals']
  
  K_all = K_diploid(simulation["Ks"]["K_all"], maternals, paternals)
  K_cas = K_diploid(simulation["Ks"]["K_cas"], maternals, paternals)
  K_obs = K_diploid(simulation["Ks"]["K_obs"], maternals, paternals)
  Km = K_diploid(simulation["Ks"]["Km"], maternals, paternals)
  Km_relate = K_diploid(simulation["Ks"]["Km_relate"], maternals, paternals)
  Km_tsinfer = K_diploid(simulation["Ks"]["Km_tsinfer"], maternals, paternals)
  Km_tsdate = K_diploid(simulation["Ks"]["Km_tsdate"], maternals, paternals)
  
  a = []
  b = []
  c = []
  d = []
  e = []
  f = []
  g = []
  for i in range(repeats):
    #print(i)
    tests = np.random.choice(int(N/2), math.floor(N * 0.125), replace = False)
    tests.sort()
    trains = [i for i in range(int(N/2)) if i not in tests]
    y_train = y_diploid[trains]
    y_test = y_diploid[tests]
    
    y_ = BLUP(K_all, y_train, trains, tests, h2 = 0.9)
    a.append(np.corrcoef(y_, y_test)[0, 1])
    y_ = BLUP(K_cas, y_train, trains, tests, h2 = 0.9)
    b.append(np.corrcoef(y_, y_test)[0, 1])
    y_ = BLUP(K_obs, y_train, trains, tests, h2 = 0.9)
    c.append(np.corrcoef(y_, y_test)[0, 1])
    y_ = BLUP(Km, y_train, trains, tests, h2 = 0.9)
    d.append(np.corrcoef(y_, y_test)[0, 1])
    y_ = BLUP(Km_relate, y_train, trains, tests, h2 = 0.9)
    e.append(np.corrcoef(y_, y_test)[0, 1])
    y_ = BLUP(Km_tsinfer, y_train, trains, tests, h2 = 0.9)
    f.append(np.corrcoef(y_, y_test)[0, 1])
    y_ = BLUP(Km_tsdate, y_train, trains, tests, h2 = 0.9)
    g.append(np.corrcoef(y_, y_test)[0, 1])
  
  a = np.array(a)
  b = np.array(b)
  c = np.array(c)
  d = np.array(d)
  e = np.array(e)
  f = np.array(f)
  g = np.array(g)
  blup_diploid = {"K_all":a, "K_cas":b, "K_obs":c, "Km":d, "Km_relate":e, "Km_tsinfer":f, "Km_tsdate":g}
  simulation["tests"]['blup_diploid'] = blup_diploid
'''

def h_estimate(K, y):
  N = y.shape[0]
  y_norm = y - y.mean(); y_norm = y_norm / y_norm.std()
  K2 = np.diag(np.dot(K, K)).sum()
  yTKy = np.dot(np.dot(np.transpose(y_norm), K), y_norm)
  yTy = np.dot(np.transpose(y_norm), y_norm)
  h2g = (yTKy - yTy) / (K2 - N)
  
  I = np.identity(N)
  Sigma = h2g * K + (1 - h2g) * I
  tmp = np.dot(Sigma, (K - I))
  tmp = np.diag(np.dot(tmp, tmp)).sum()
  variance = np.sqrt(tmp) / ((K2 - N) * h2g)
  return h2g, variance #needs change its usage accordingly!!!

def test(simulation, repeats = 100):
  N = simulation["parameters"]["N"]
  y = simulation["phenotypes"]["y"]
  Ks = getKs(simulation)
  
  diags = np.diag_indices(N)
  non_diags = np.where(~np.eye(N,dtype=bool))
  
  table = {}
  for key in Ks.keys():
    table[key] = Ks[key][non_diags].flatten()
  
  table = pd.DataFrame(data=table)
  corr = table.corr(method ='pearson')
  
  h_estimation = {}
  for key in Ks.keys():
    h_estimation[key] = h_estimate(Ks[key], y)
  
  #p_imputation = phenotype_impute(simulation, repeats)
  p_imputation = {}
  
  simulation["tests"] = {"corr":corr, 'h_estimation':h_estimation, 'p_imputation':p_imputation}


def summary(simulation):
  summ = "==========\nparameters \n==========\n"
  summ += "\n".join([str(x) + "\t" + str(simulation["parameters"][x]) for x in simulation["parameters"].keys()]) + "\n"
  
  summ += "==========\nK matrix correlations \n==========\n"
  summ += str(simulation["tests"]["corr"]) + "\n"
  
  summ += "==========\nheritability estimation \n==========\n"
  for key in simulation["tests"]["h_estimation"].keys():
    tmp = simulation["tests"]["h_estimation"][key]
    summ += key + "\t" + str(round(tmp[0], 4)) + " +- " + str(round(tmp[1], 4)) + "\n"
  
  summ += "==========\nBLUP accuracy \n==========\n"
  for key in simulation["tests"]["p_imputation"].keys():
    tmp = simulation["tests"]["p_imputation"][key]
    summ += key + "\t" + str(round(tmp.mean(), 4)) + " +- " + str(round(tmp.std(), 4)) + "\n"
  
  return(summ)

