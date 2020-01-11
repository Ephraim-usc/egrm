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

def test(simulation, repeats = 1000):
  N = simulation["parameters"]["N"]
  K_cas = simulation["Ks"]["K_cas"]
  K_obs = simulation["Ks"]["K_obs"]
  Km = simulation["Ks"]["Km"]
  Km_relate = simulation["Ks"]["Km_relate"]
  
  diags = np.diag_indices(N)
  non_diags = np.where(~np.eye(N,dtype=bool))
  
  table = {"K_cas":K_cas[non_diags].flatten(), "K_obs":K_obs[non_diags].flatten(),
           "Km":Km[non_diags].flatten(), "Km_relate":Km_relate[non_diags].flatten()}
  
  table = pd.DataFrame(data=table)
  corr = table.corr(method ='pearson')
  
  y = simulation["phenotypes"]["y"]
  a = []
  b = []
  c = []
  d = []
  for i in range(1000):
    tests = np.random.choice(N, math.floor(N * 0.25), replace = False)
    tests.sort()
    trains = [i for i in range(N) if i not in tests]
    y_train = y[trains]
    y_test = y[tests]
    
    y_ = BLUP(K_cas, y_train, trains, tests, h2 = 0.9)
    a.append(np.corrcoef(y_, y_test)[0, 1])
    y_ = BLUP(K_obs, y_train, trains, tests, h2 = 0.9)
    b.append(np.corrcoef(y_, y_test)[0, 1])
    y_ = BLUP(Km, y_train, trains, tests, h2 = 0.9)
    c.append(np.corrcoef(y_, y_test)[0, 1])
    y_ = BLUP(Km_relate, y_train, trains, tests, h2 = 0.9)
    d.append(np.corrcoef(y_, y_test)[0, 1])
  
  a = np.array(a)
  b = np.array(b)
  c = np.array(c)
  d = np.array(d)
  blup = {"K_cas":a, "K_obs":b, "Km":c, "Km_relate":d}
  simulation["tests"] = {"corr":corr, "blup":blup}



