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

### run workflow
out = args["out"]
os.system("mkdir " + out)
os.chdir(out)
del args['out']

relate = args["relate"]
del args['relate']

simulation = simulate(**args)
make_diploid(simulation)

M = simulation["hapdata"]["M"]
cass = simulation["phenotypes"]["cass"]
obss = simulation["observations"]["obss"]
write_plink(simulation, obss, "observed")
write_plink(simulation, cass, "causal")

gcta64("observed")
gcta64("causal")

if relate:
  write_relate("observed")
  relate("observed")
