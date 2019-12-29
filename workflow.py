### loading packages
import argparse
import os

### loading files
exec(open("crm/crm.py").read())
exec(open("crm/read_and_write.py").read())
exec(open("crm/simulation.py").read())
exec(open("crm/blup.py").read())

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

def printf(string):
  out_file = open(name, "a")
  out_file.write(str(string) + "\n")
  out_file.close()

### run workflow
printf("simulating ...")
simulation = simulate(**args)
make_diploid(simulation)

trees = simulation["hapdata"]["tree_sequence"]
N = simulation["hapdata"]["N"]
M = simulation["hapdata"]["M"]
M_cas = simulation["phenotypes"]["M_cas"]
M_obs = simulation["observations"]["M_obs"]
M_5 = simulation["observations"]["M_5"]
cass = simulation["phenotypes"]["cass"]
obss = simulation["observations"]["obss"]

printf("computing Km ...")
flags_obs = get_flags(trees, obss)
Km = getEK_trees(trees, flags_obs)

if run_gcta:
  printf("running gcta ...")
  os.mkdir(name + "_gcta")
  os.chdir(name + "_gcta")
  write_plink(simulation, obss, "observed")
  write_plink(simulation, cass, "causal")
  gcta64("observed")
  gcta64("causal")
  os.chdir("..")

if run_relate:
  printf("running relate ...")
  os.mkdir(name + "_relate")
  os.chdir(name + "_relate")
  write_relate(simulation, obss, "observed")
  relate("observed")
  trees_relate = tskit.load("observed.trees")
  printf("computing Km_relate ...")
  Km_relate = getEK_trees(trees_relate)
  os.chdir("..")

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
  
  if run_relate == False:
    continue
  
  y_ = BLUP(Km_relate, y_train, trains, tests, h2 = 0.9)
  d.append(np.corrcoef(y_, y_test)[0, 1])

### output
a = np.array(a)
b = np.array(b)
c = np.array(c)
d = np.array(d)

printf(args)
printf("M_cas / M = " + str(np.array(M_cas/M).round(2)))
printf("M_obs / M_5 = " + str(np.array(M_obs/M_5).round(2)))

printf("true number of trees = " + str(trees.num_trees))
printf("relate inferred number of trees = " + str(trees_relate.num_trees))

diags = np.diag_indices(N)
non_diags = np.where(~np.eye(N,dtype=bool))

table = {"K_cas":K_cas[non_diags].flatten(), "K_obs":K_obs[non_diags].flatten(),
         "Km":Km[non_diags].flatten(), "Km_relate":Km_relate[non_diags].flatten()}

table = pd.DataFrame(data=table)
printf(table.corr(method ='pearson'))

printf("K_cas BLUP: " + str(a.mean().round(3)) + " +- " + str(a.std().round(3)))
printf("K_obs BLUP: " + str(b.mean().round(3)) + " +- " + str(b.std().round(3)))
printf("Km BLUP: " + str(c.mean().round(3)) + " +- " + str(c.std().round(3)))
printf("Km_relate BLUP: " + str(d.mean().round(3)) + " +- " + str(d.std().round(3)))

