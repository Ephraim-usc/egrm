### loading packages
import argparse
import os

### loading files
exec(open("crm/crm.py").read())
exec(open("crm/read_and_write.py").read())
exec(open("crm/simulation.py").read())

### parse arguments
parser=argparse.ArgumentParser()

parser.add_argument('--out', type = str, help='output directory')
parser.add_argument('--name', type = str, help='output file name')
parser.add_argument('--gcta', help='run gcta analysis', action='store_true')
parser.add_argument('--relate', help='run relate tree reconstruction', action='store_true')
parser.add_argument('--l', type = int, help='chromosome length')
parser.add_argument('--N', type = int, help='population size')
parser.add_argument('--mutation_rate', type = float, help='mutation rate')
parser.add_argument('--recomb_rate', type = float, help='recombination rate')
parser.add_argument('--cas_ratio', type = float, help='M_cas / M')
parser.add_argument('--obs_ratio', type = float, help='M_obs / M_5')
parser.add_argument('--h2g', type = float, help='heritability')
parser.add_argument('--Alpha', type = float, help='Alpha -- causal probability parameter')
parser.add_argument('--Beta', type = float, help='Beta -- observation probability parameter')

args = vars(parser.parse_args())
args = dict((k,v) for k,v in args.items() if v is not None)

out = args["out"]
os.chdir(out)
del args['out']

name = args["name"]
del args['name']

run_gcta = args["gcta"]
del args['gcta']

run_relate = args["relate"]
del args['relate']

### run workflow
simulation = simulate(**args)
make_diploid(simulation)

M = simulation["hapdata"]["M"]
cass = simulation["phenotypes"]["cass"]
obss = simulation["observations"]["obss"]

if run_gcta:
  write_plink(simulation, obss, name + "observed")
  write_plink(simulation, cass, name + "causal")
  gcta64(name + "observed")
  gcta64(name + "causal")

if run_relate:
  write_relate(simulation, obss, name + "observed_relate")
  relate(name + "observed_relate")
  trees_relate = tskit.load(name + "observed_relate.trees")
  K_relate = getEK_trees(trees_relate)

### run BLUP
Z_cas = simulation["phenotypes"]["Z_cas"]
Z_obs = simulation["observations"]["Z_obs"]
K_cas = np.dot(Z_cas, np.transpose(Z_cas)) / Z_cas.shape[1]
K_obs = np.dot(Z_obs, np.transpose(Z_obs)) / Z_obs.shape[1]
N = Z_cas.shape[0]
y = simulation["phenotypes"]["y"]

a = []
b = []
c = []
for i in range(100):
  tests = np.random.choice(N, math.floor(N * 0.25), replace = False)
  tests.sort()
  trains = [i for i in range(N) if i not in tests]
  y_train = y[trains]
  y_test = y[tests]
  
  y_ = BLUP(K_cas, y_train, trains, tests, h2 = 0.9)
  a.append(np.corrcoef(y_, y_test)[0, 1])
  
  y_ = BLUP(K_obs, y_train, trains, tests, h2 = 0.9)
  b.append(np.corrcoef(y_, y_test)[0, 1])
  
  y_ = BLUP(K_relate, y_train, trains, tests, h2 = 0.9)
  c.append(np.corrcoef(y_, y_test)[0, 1])

a = np.array(a)
b = np.array(b)
c = np.array(c)

