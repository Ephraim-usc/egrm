### loading packages
import argparse

### loading files
exec(open("crm/crm.py").read())
exec(open("crm/read_and_write.py").read())
exec(open("crm/simulation.py").read())

###
parser=argparse.ArgumentParser()

parser.add_argument('--Alpha', help='Alpha -- causal probability parameter')
parser.add_argument('--Beta', help='Beta -- observation probability parameter')

args=parser.parse_args()

print(args)
print(args["Alpha"])
