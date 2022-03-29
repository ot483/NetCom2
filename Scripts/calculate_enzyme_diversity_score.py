#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 14:45:08 2022

@author: ofir
"""




import sys
import pandas as pd
#import re
#import numpy as np
#from multiprocessing import Pool
#from ete3 import NCBITaxa
import skbio
from skbio.diversity.alpha import simpson


output = sys.argv[1]
BaseFolder = sys.argv[2]


#R script
data = pd.read_csv(BaseFolder + output, sep = '\t')
dominance_cutoff = 0.4

cols = list(data.columns)
cols = [i for i in cols if not i in ["EC", "taxa"]]

data[cols] = data[cols].astype(int)

Unique_EC_List = data["EC"].unique()

All_simpson = []
#All_shannon = []
Shannon_temp = []
Simpson_temp = []
Sum_table = []
All_sumTables = {}

for j in Unique_EC_List:
    df_temp = data[data["EC"] == j][['taxa']+cols]
    df_temp = df_temp.set_index("taxa")
    Simpson_temp.append(j)
    
    for i in cols:       
        SimpsonScore = simpson(df_temp[i].values.tolist())
  
        #Most frequent organism
        mostFreqOrg = df_temp[i].idxmax()
        Dominance = 1-SimpsonScore
        Simpson_temp.append(SimpsonScore)            
        Simpson_temp.append(Dominance)
        Simpson_temp.append(mostFreqOrg)
             
    All_simpson.append(Simpson_temp)
    Simpson_temp = []
    
    print("EC "+j+" done")

newCols = []
newCols.append("Enzyme")
DominanceCols = []
SimpsonCols = []
MostFreqCols = []

for i in cols:
    newCols.append("SimpsonScore_"+i)
    newCols.append("Dominance_"+i)
    newCols.append("mostFreqOrg_"+i)
    MostFreqCols.append("mostFreqOrg_"+i)
    SimpsonCols.append("SimpsonScore_"+i)
    DominanceCols.append("Dominance_"+i)

final_df = pd.DataFrame(data=All_simpson, columns=newCols)

allEnzymesList = final_df["Enzyme"].values.tolist()
allEnzymesList = list(set(allEnzymesList))

final_df["DominanceMean"] = final_df[DominanceCols].mean(axis=1)

final_df.to_csv(BaseFolder+"RscriptOutput.csv")

final_df = final_df[final_df["DominanceMean"] >= dominance_cutoff]

final_df.to_csv(BaseFolder+"RscriptOutput_filter_dominance.csv")

final_df["AllIdent"] = final_df[MostFreqCols].all(axis=1)

final_df = final_df[final_df["AllIdent"] == True]

final_df = final_df[final_df[newCols[-1]] != "Unidentified"][newCols]

final_df.to_csv(BaseFolder+"RscriptOutput_filter_IdentMostFreq.csv")

final_df_byEnzyme = final_df[["Enzyme", MostFreqCols[-1]]].groupby(MostFreqCols[-1])["Enzyme"].apply(list).to_frame()
final_df_byEnzyme = final_df_byEnzyme.reset_index()
final_df_byEnzyme.columns = ["Organism", "Enzyme"]

final_df_byEnzyme.to_csv(BaseFolder+"RscriptOutput_filter_IdentMostFreqbyEnzyme.csv")

#Add knockout
final_df_byEnzyme_knockout = final_df_byEnzyme.copy()

allEnzs = final_df_byEnzyme_knockout.values.tolist()

for i, vali in enumerate(allEnzs):
    allEnzs[i][1] = [j for j in allEnzymesList if not j in vali[1] ]

final_df_byEnzyme_knockout = pd.DataFrame(data=allEnzs, columns=final_df_byEnzyme.columns)

def fix_output(x):
    return " ".join(x).strip(" ")

final_df_byEnzyme.loc[-1] = ['all', allEnzymesList]  # adding a row
final_df_byEnzyme.index = final_df_byEnzyme.index + 1  # shifting index
final_df_byEnzyme.sort_index(inplace=True) 

final_df_byEnzyme["Enzyme"] = final_df_byEnzyme["Enzyme"].apply(fix_output)
final_df_byEnzyme.to_csv(BaseFolder+"ECs_with_dominant_taxon.txt", sep=" ", index=False, header=False)


final_df_byEnzyme_knockout.loc[-1] = ['all', allEnzymesList]  # adding a row
final_df_byEnzyme_knockout.index = final_df_byEnzyme_knockout.index + 1  # shifting index
final_df_byEnzyme_knockout.sort_index(inplace=True) 

final_df_byEnzyme_knockout["Enzyme"] = final_df_byEnzyme_knockout["Enzyme"].apply(fix_output)
final_df_byEnzyme_knockout.to_csv(BaseFolder+"ECs_with_dominant_taxon_knockout.txt", sep=" ", index=False, header=False)

lines = []
infile = open(BaseFolder+"ECs_with_dominant_taxon.txt", "r")
for line in infile.readlines():
    line = str(line).replace("\"", "") # remove the newline at the end
    #print(line) # or do something else
    lines.append(line)
#write output files

open(BaseFolder+"ECs_with_dominant_taxon.txt", 'w').writelines(lines)

lines = []
infile = open(BaseFolder+"ECs_with_dominant_taxon_knockout.txt", "r")
for line in infile.readlines():
    line = str(line).replace("\"", "") # remove the newline at the end
    #print(line) # or do something else
    lines.append(line)
#write output files

open(BaseFolder+"ECs_with_dominant_taxon_knockout.txt", 'w').writelines(lines)
####




