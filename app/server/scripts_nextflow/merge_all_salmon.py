import os
import pandas as pd 
import sys


array = [sys.argv[i] for i in range(1,len(sys.argv)-2)]
new_data_path1 = sys.argv[-2]
new_data_path2 = sys.argv[-1]
#print(array)

path = ""
bool_existence = False
for i,val in enumerate(array):
    path = val + "salmon/quant.sf"
    if os.path.exists(path):
        template_df1 = pd.read_csv(path, sep = "\t",header = 0)
        template_df1 = template_df1.iloc[:,0:3]
        template_df2 = pd.read_csv(path, sep = "\t",header = 0)
        template_df2 = template_df2.iloc[:,0:3]
        bool_existence = True
        break
    else:
        execute = False
        #print("None of the quant.sf in this path exist")

if bool_existence:
    for i,val in enumerate(array):
        path = val + "salmon/quant.sf"  
        if os.path.exists(path):  
            folder_name = val.split("/")
            folder_name = folder_name[-2]
            print("Folder name 2: ",folder_name)
            df_i = pd.read_csv(path, sep="\t", header = 0)
            colnames = list(df_i.columns)
            colnames[-1] = folder_name
            colnames[-2] = folder_name
            print("Columns 2: ",colnames)
            df_i.columns = colnames
            template_df1 = pd.concat([template_df1,df_i.iloc[:,-2]], axis = 1)
            template_df2 = pd.concat([template_df2,df_i.iloc[:,-1]], axis = 1)
        else:
            folder_name = val.split("/")
            folder_name = folder_name[-2]
            zero_list1 = [0 for i in range(len(template_df1.iloc[:,1]))]
            zero_list2 = [0 for i in range(len(template_df2.iloc[:,1]))]
            template_df1 = pd.concat([template_df1,pd.DataFrame(zero_list1)], axis = 1)
            template_df2 = pd.concat([template_df2,pd.DataFrame(zero_list2)], axis = 1)
            colnames1 = list(template_df1.columns)
            colnames1[-1] = folder_name
            template_df1.columns = colnames1
            colnames2 = list(template_df2.columns)
            colnames2[-1] = folder_name
            template_df2.columns = colnames2
    template_df1["Name"] = [i.split("|")[0] for i in list(template_df1["Name"])]
    template_df2["Name"] = [i.split("|")[0] for i in list(template_df2["Name"])]
    template_df1.to_csv(new_data_path1, sep="\t", index = True)
    template_df2.to_csv(new_data_path2, sep="\t", index = True)
        
     
