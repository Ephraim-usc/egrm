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

### parameters
out = "tmp"
log_file = "tmp.log"
l = 32000000
N = 1000
mutation_rate = 1e-8
recomb_rate = 1e-8
cas_ratio = 0.1
obs_ratio = 0.2
Alpha = -0.5
Beta = 1
h2g = 1.0

### out-of-Africa parameters
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



    

    
