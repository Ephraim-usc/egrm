### loading packages
import argparse

### loading files
exec(open("crm/crm.py").read())
exec(open("crm/read_and_write.py").read())
exec(open("crm/simulation.py").read())

###
parser=argparse.ArgumentParser()

parser.add_argument('--l', type = int, help='chromosome length')
parser.add_argument('--N', type = int, help='population size')
parser.add_argument('--mutation_rate', type = float, help='mutation rate')
parser.add_argument('--recomb_rate', type = float, help='recombination rate')
parser.add_argument('--cas_ratio', type = float, help='M_cas / M')
parser.add_argument('--obs_ratio', type = float, help='M_obs / M_5')
parser.add_argument('--h2g', type = float, help='heritability')
parser.add_argument('--Alpha', type = float, help='Alpha -- causal probability parameter')
parser.add_argument('--Beta', type = float, help='Beta -- observation probability parameter')

args=parser.parse_args()

simulation = simulate(**vars(args))

print(simulation["hapdata"]["genotype_martix"].shape)
