import os
import pandas as pd 
import sys


array = [sys.argv[i] for i in range(1,len(sys.argv)-1)]
new_data_path = sys.argv[-1]
print(array)

path = ""
bool_existence = False
for i,val in enumerate(array):
    path = val + "merged_fc/merged_fc.csv"
    if os.path.exists(path):
        template_df = pd.read_csv(path, sep = "\t",header = 1)
        template_df = template_df.iloc[:,0]
        bool_existence = True
        break
    else:
        execute = False
        print("None of the .csv in this path exist")

if bool_existence:
    for i,val in enumerate(array):
        path = val + "merged_fc/merged_fc.csv"  
        if os.path.exists(path):  
            folder_name = val.split("/")
            folder_name = folder_name[-2]
            print("Folder name 2: ",folder_name)
            df_i = pd.read_csv(path, sep="\t", header = 1)
            colnames = list(df_i.columns)
            colnames[-1] = folder_name
            print("Columns 2: ",colnames)
            df_i.columns = colnames
            template_df = pd.concat([template_df,df_i.iloc[:,-1]], axis = 1)
        else:
            folder_name = val.split("/")
            folder_name = folder_name[-2]
            print("Folder name not existent: ", folder_name)
            zero_list = [0 for i in range(len(template_df.iloc[:]))]
            template_df = pd.concat([template_df,pd.DataFrame(zero_list)], axis = 1)
            colnames = list(template_df.columns)
            colnames[-1] = folder_name
            template_df.columns = colnames
    #print(template_df.head())

    template_df.to_csv(new_data_path, sep="\t", index = False)
        
     
       
        
        
        
        
    
    
