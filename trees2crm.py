### loading packages
import argparse
import tskit
import pandas as pd

### loading files
exec(open("crm/crm.py").read())
exec(open("crm/read_and_write.py").read())

### parse arguments
parser=argparse.ArgumentParser()
parser.add_argument('--input', type = str, help='input and output file name')

args = vars(parser.parse_args())
args = dict((k,v) for k,v in args.items() if v is not None)

input = args["input"]
output = input[:-6]

### 
trees = tskit.load(input)
crm = getEK_trees(trees)

write_crm(crm, output)
