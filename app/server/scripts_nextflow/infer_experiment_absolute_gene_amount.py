import argparse
import pandas as pd
import numpy as np
import os



opt_parser = argparse.ArgumentParser()

opt_parser.add_argument("-s", "--sample_file", dest="sample", help="Insert a sample file to add names to", metavar="FILE")
opt_parser.add_argument("-m", "--metadata_file",dest="metadata", help="Insert a metadata file to extract metdata from", metavar="FILE")
opt_parser.add_argument("-o", "--output_path",dest="output", help="Insert a template file to extract names from", metavar="FILE")

options = opt_parser.parse_args()

sample = options.sample
metadata = options.metadata
output_path = options.output

sample_df = pd.read_csv(sample,header = 0, index_col = 0, sep = "\t")
#print(sample_df)
if not os.path.exists(output_path):
    output_df = pd.DataFrame()
else:
    output_df = pd.read_csv(output_path,header=0, sep = "\t")

samplenames = []
genes_counted = []
for i in range(len(sample_df.iloc[0,:])):
    column = sample_df.iloc[:,i]
    #print(column[0])
    if not type(column[0]) == type(""):
        name = column.name
        samplenames.append(column.name)
        column = pd.DataFrame(column)
        #print(column[name])
        length = len(column.loc[column[name] > 0,:])
        genes_counted.append(length)

#print(samplenames)
#print(genes_counted)

input_list = [genes_counted]
new_row_df = pd.DataFrame(np.array(input_list), columns = samplenames)

#print(new_row_df)

output_df = pd.concat([output_df,new_row_df])
output_df = output_df.reset_index(drop = True)
#print(output_df)
output_df.to_csv(output_path, sep="\t", index = 0)
#print(output_df)








    


