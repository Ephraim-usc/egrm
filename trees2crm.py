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
output = input[:-6] + ".grm"

### 
trees = tskit.load("input")
crm = getEK_trees(trees)

crm_df = pd.DataFrame(data = crm)
crm_df["row"] = range(crm.shape[0])
crm_melted = pd.melt(crm_df, id_vars=['row'], var_name='column', value_name='crm')


crm_melted = pd.melt(crm_df)
