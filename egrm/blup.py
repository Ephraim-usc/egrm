import numpy as np
import pandas as pd
import math

def BLUP(K, y_train, trains, tests, h2 = 0.9):
  N_train = len(trains)
  
  I = (1/h2 - h2) * np.identity(N_train)
  V = K[trains, trains] + I
  
  Kt = K[tests, :][:, trains]
  
  Vi = np.linalg.inv(V)
  y_ = np.dot(Kt, np.dot(Vi, y_train))
  return y_

def phenotype_impute(simulation, repeats = 100):
  y = simulation['phenotypes']['y']
  N = simulation["parameters"]["N"]
  K_all = simulation["Ks"]["K_all"]
  K_cas = simulation["Ks"]["K_cas"]
  K_obs = simulation["Ks"]["K_obs"]
  Km = simulation["Ks"]["Km"]
  Km_relate = simulation["Ks"]["Km_relate"]
  Km_tsinfer = simulation["Ks"]["Km_tsinfer"]
  
  a = []
  b = []
  c = []
  d = []
  e = []
  f = []
  for i in range(repeats):
    #print(i)
    tests = np.random.choice(N, math.floor(N * 0.25), replace = False)
    tests.sort()
    trains = [i for i in range(N) if i not in tests]
    y_train = y[trains]
    y_test = y[tests]
    
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
  
  a = np.array(a)
  b = np.array(b)
  c = np.array(c)
  d = np.array(d)
  e = np.array(e)
  f = np.array(f)
  blup = {"K_all":a, "K_cas":b, "K_obs":c, "Km":d, "Km_relate":e, "Km_tsinfer":f}
  simulation["tests"]['blup'] = blup

def K_diploid(K, maternals, paternals):
  K1 = K[maternals, :][:, maternals]
  K2 = K[maternals, :][:, paternals]
  K3 = K[paternals, :][:, maternals]
  K4 = K[paternals, :][:, paternals]
  return (K1 + K2 + K3 + K4)/2

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
  
  a = []
  b = []
  c = []
  d = []
  e = []
  f = []
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
  
  a = np.array(a)
  b = np.array(b)
  c = np.array(c)
  d = np.array(d)
  e = np.array(e)
  f = np.array(f)
  blup_diploid = {"K_all":a, "K_cas":b, "K_obs":c, "Km":d, "Km_relate":e, "Km_tsinfer":f}
  simulation["tests"]['blup_diploid'] = blup_diploid

def h_estimate(K, y, N):
  y_norm = y - y.mean(); y_norm = y_norm / y_norm.std()
  K2 = np.diag(np.dot(K, K)).sum()
  yTKy = np.dot(np.dot(np.transpose(y_norm), K), y_norm)
  yTy = np.dot(np.transpose(y_norm), y_norm)
  h2g = (yTKy - yTy) / (K2 - N)
  
  I = np.identity(N)
  Sigma = h2g * K + (1 - h2g) * I
  tmp = np.dot(Sigma, (K - I))
  tmp = np.diag(np.dot(tmp, tmp)).sum()
  variance = tmp / ((K2 - N) * h2g)
  return h2g, variance #needs change its usage accordingly!!!

def test(simulation, repeats = 1000):
  N = simulation["parameters"]["N"]
  K_all = simulation["Ks"]["K_all"]
  K_cas = simulation["Ks"]["K_cas"]
  K_obs = simulation["Ks"]["K_obs"]
  Km = simulation["Ks"]["Km"]
  Km_relate = simulation["Ks"]["Km_relate"]
  Km_tsinfer = simulation["Ks"]["Km_tsinfer"]
  
  diags = np.diag_indices(N)
  non_diags = np.where(~np.eye(N,dtype=bool))
  
  table = {"K_all":K_all[non_diags].flatten(),
           "K_cas":K_cas[non_diags].flatten(), 
           "K_obs":K_obs[non_diags].flatten(),
           "Km":Km[non_diags].flatten(), 
           "Km_relate":Km_relate[non_diags].flatten(),
           "Km_tsinfer":Km_tsinfer[non_diags].flatten()}
  
  table = pd.DataFrame(data=table)
  corr = table.corr(method ='pearson')
  
  y = simulation["phenotypes"]["y"]
  h_estimation = {'K_all':h_estimate(K_all, y, N), 
                  'K_cas':h_estimate(K_cas, y, N), 
                  'K_obs':h_estimate(K_obs, y, N),
                  'Km':h_estimate(Km, y, N),
                  'Km_relate':h_estimate(Km_relate, y, N),
                  'Km_tsinfer':h_estimate(Km_tsinfer, y, N)} 
  
  simulation["tests"] = {"corr":corr, 'h_estimation':h_estimation}
  phenotype_impute(simulation, repeats)
  phenotype_impute_diploid(simulation, repeats)


def summary(simulation):
  summ = "==========\nparameters \n==========\n"
  summ += "\n".join([str(x) + "\t" + str(simulation["parameters"][x]) for x in simulation["parameters"].keys()]) + "\n"
  
  summ += "==========\nK matrix correlations \n==========\n"
  summ += str(simulation["tests"]["corr"]) + "\n"
  
  summ += "==========\nheritability estimation \n==========\n"
  summ += "K_all\t" + str(round(simulation["tests"]["h_estimation"]["K_all"], 4)) + "\n"
  summ += "K_cas\t" + str(round(simulation["tests"]["h_estimation"]["K_cas"], 4)) + "\n"
  summ += "K_obs\t" + str(round(simulation["tests"]["h_estimation"]["K_obs"], 4)) + "\n"
  summ += "Km\t" + str(round(simulation["tests"]["h_estimation"]["Km"], 4)) + "\n"
  summ += "Km_relate\t" + str(round(simulation["tests"]["h_estimation"]["Km_relate"], 4)) + "\n"
  summ += "Km_tsinfer\t" + str(round(simulation["tests"]["h_estimation"]["Km_tsinfer"], 4)) + "\n"
  
  summ += "==========\nBLUP accuracy \n==========\n"
  tmp = simulation["tests"]["blup"]["K_all"]
  summ += "K_all\t" + str(round(tmp.mean(), 4)) + " +- " + str(round(tmp.std(), 4)) + "\n"
  tmp = simulation["tests"]["blup"]["K_cas"]
  summ += "K_cas\t" + str(round(tmp.mean(), 4)) + " +- " + str(round(tmp.std(), 4)) + "\n"
  tmp = simulation["tests"]["blup"]["K_obs"]
  summ += "K_obs\t" + str(round(tmp.mean(), 4)) + " +- " + str(round(tmp.std(), 4)) + "\n"
  tmp = simulation["tests"]["blup"]["Km"]
  summ += "Km\t" + str(round(tmp.mean(), 4)) + " +- " + str(round(tmp.std(), 4)) + "\n"
  tmp = simulation["tests"]["blup"]["Km_relate"]
  summ += "Km_relate\t" + str(round(tmp.mean(), 4)) + " +- " + str(round(tmp.std(), 4)) + "\n"
  tmp = simulation["tests"]["blup"]["Km_tsinfer"]
  summ += "Km_tsinfer\t" + str(round(tmp.mean(), 4)) + " +- " + str(round(tmp.std(), 4)) + "\n"
  
  summ += "==========\ndiploid BLUP accuracy \n==========\n"
  tmp = simulation["tests"]["blup_diploid"]["K_all"]
  summ += "K_all\t" + str(round(tmp.mean(), 4)) + " +- " + str(round(tmp.std(), 4)) + "\n"
  tmp = simulation["tests"]["blup_diploid"]["K_cas"]
  summ += "K_cas\t" + str(round(tmp.mean(), 4)) + " +- " + str(round(tmp.std(), 4)) + "\n"
  tmp = simulation["tests"]["blup_diploid"]["K_obs"]
  summ += "K_obs\t" + str(round(tmp.mean(), 4)) + " +- " + str(round(tmp.std(), 4)) + "\n"
  tmp = simulation["tests"]["blup_diploid"]["Km"]
  summ += "Km\t" + str(round(tmp.mean(), 4)) + " +- " + str(round(tmp.std(), 4)) + "\n"
  tmp = simulation["tests"]["blup_diploid"]["Km_relate"]
  summ += "Km_relate\t" + str(round(tmp.mean(), 4)) + " +- " + str(round(tmp.std(), 4)) + "\n"
  tmp = simulation["tests"]["blup_diploid"]["Km_tsinfer"]
  summ += "Km_tsinfer\t" + str(round(tmp.mean(), 4)) + " +- " + str(round(tmp.std(), 4)) + "\n"
  
  return(summ)

