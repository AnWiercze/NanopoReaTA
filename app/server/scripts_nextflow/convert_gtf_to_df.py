import argparse
import pandas as pd
import numpy as np 
import gtfparse as gtf_parse

opt_parser = argparse.ArgumentParser()
opt_parser.add_argument("-i", "--input", dest="input_file", help="Insert a gtf file to parse", metavar="FILE")
opt_parser.add_argument("-o", "--output",dest="output_file", help="Insert a gpath for the output file", metavar="FILE")
options = opt_parser.parse_args()



input_file = options.input_file
#print(input_file)
df = gtf_parse.read_gtf(input_file)

#print(df.head(20))
#print(type(df))
#print(len(df))

output_file = options.output_file
#print(output_file)
pd.DataFrame.to_csv(df,output_file)
