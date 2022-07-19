#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 20 08:49:55 2022

@author: ofir
"""


import sys, os
import pandas as pd
import re
import numpy as np
from multiprocessing import Pool
from ete3 import NCBITaxa
import skbio
from skbio.diversity.alpha import simpson



def create_temp_file(input_, tempfile):
    with open(input_) as f:
        with open(tempfile, 'w') as t:
            data = f.read()
            data = re.sub(r'([0-9][0-9]*_[0-9][0-9]*)_([0-9][0-9]*\t)', r'\1\t\2', data)
            data = re.sub(' ', '_', data)
            t.write(data)

def most_common(lst):
    current = max(set(lst), key=lst.count)
    maxnum = lst.count(current)
    listmax = []
    for el in set(lst):
        if lst.count(el) == maxnum:
            listmax.append(el)
    if len(listmax) == 1 :
        return listmax[0]
    else:
        return 'even'

def manipulate(cont):
    results_list = []
    #global number_current_contig
    #number_current_contig += 1
    #print(str(number_current_contig) + " / " + str(contnum) )
    #print(("Starting to work on contig {} of {}").format(str(), str(contnum)))
    contdata = data.where(data['contig'] == cont).dropna()
    smalldata = contdata['annotation'].str.split(';', expand=True)
    smalldata.replace('', np.nan, inplace=True)
    smalldata.replace('Not_assigned', np.nan, inplace=True)
    most_freq_tax = ''
    common_tax = ''
    levelnum = len(smalldata.columns)
     #I set a parameter (diff) that switch to 1 from 0 if the intra-contig taxonomy differs, so the program stop adding things to the common taxonomy column.
    diff = 0
    for num in range(levelnum):
        anlist = smalldata[num].dropna().tolist()
        if len(set(anlist)) == 1:                        
            most_freq_tax += (';' + anlist[0])
            if diff == 0:
                common_tax += (';' + anlist[0])
        elif len(set(anlist)) > 1:
            toadd = most_common(anlist)
            if toadd != 'even':
                most_freq_tax += (';' + toadd)
                smalldata = smalldata[smalldata[num] == toadd]
                diff = 1
            else:
                break
    if most_freq_tax == common_tax:
        #o.write(cont + '\t' + most_freq_tax + '\t' + common_tax + '\t' + 'Yes' + '\n')
        results_list.append([cont + '\t' + most_freq_tax + '\t' + common_tax + '\t' + 'Yes' + '\n'])
    else:
        #o.write(cont + '\t' + most_freq_tax + '\t' + common_tax + '\t' + 'No' + '\n')
        results_list.append([cont + '\t' + most_freq_tax + '\t' + common_tax + '\t' + 'No' + '\n'])
    
    return results_list

def get_last_taxonomy(string):
    x = string.split(";")
    return x[-1]

def get_tax_rank(string):
    x = string.replace('_', ' ').split(", ")
    for i in range(len(x)):
        try:
            taxid = ncbi.get_name_translator([x[i]])[x[i]]
            rank = ncbi.get_rank(taxid)[taxid[0]]
            if rank != 'no rank' :
                break
        except KeyError:
            continue
    try:
        rank
    except NameError:
        rank = "Unknown_rank"
    return rank.replace(' ', '_')

def get_level_taxonomy(string, level):
    x = string.split(";")
    for i in x:
        if get_tax_rank(i) == level :
            result = i
            break
        if level != 'order':
            if get_tax_rank(i) == 'order':
                result = 'Unidentified_' + str(i)
    try:
        result
    except NameError:
        result = 'Unidentified'
    return result







input1 = sys.argv[1]
input2 = sys.argv[2]
input3 = sys.argv[3]
BaseFolder = sys.argv[4]

#input1 = "NTC_G210-ReadName_to_TaxonPath.txt"
#input2 = "NTC_G210.count_tab.matrix"
#input3 = "NTC_G210-ReadName_to_KeggName_Ultimate.txt"
#BaseFolder = "/home/ofir/Documents/Shiri/Removal_network/Input/"


#Fix ECs - converted perl script
file1 = open(BaseFolder + input3, 'r')
file1_Lines = file1.readlines()
file1.close()

file1_Lines = [i.strip('\n').split("\t") for i in file1_Lines if "EC:" in i]

for i, vali in enumerate(file1_Lines):
    start = vali[1].find('EC:')
    end = vali[1].find(']', start)
    file1_Lines[i][1] = vali[1][start+3:end]

file1_Lines = ['\t'.join(i) for i in file1_Lines]

inputUltimate_fixed = file1_Lines.copy()



####
#Script #1

file = BaseFolder + input1
tempfile = BaseFolder + "temp_1"
output = BaseFolder + "get_contig_taxonomy_ranked.txt" 


#Check input file
try:
   with open(file) as f:
       print('I will now start to work. Feel free to do something relaxing while you wait.')
except FileNotFoundError:
    print('It seems that the input file you specified does not exist. I am very sorry, but I cannot proceed :(.')
    sys.exit()

#At first I create the temporary file with the data formatted as I need
create_temp_file(file, tempfile)

data = pd.read_csv(tempfile, names=['contig', 'protein', 'annotation'], sep = '\t')
contlist = set(data['contig'].tolist())
contnum = len(contlist)

results_list_final = Pool().map(manipulate, list(contlist))

results_list_final_flatten = [ent for sublist in results_list_final for ent in sublist]
results_list_final_flatten = [['Contig' + '\t' + 'most_frequent_taxonomy' + '\t' + 'common_taxonomy' + '\t' + 'Same_taxonomy?' + '\n']]+ results_list_final_flatten

o = open(output, "w")
for element in results_list_final_flatten:    
    o.write(element[0])
o.close()



ncbi = NCBITaxa()

data = pd.read_csv(output, sep='\t')

data['taxonomy_rank'] = data['most_frequent_taxonomy'].apply(lambda x: get_tax_rank(get_last_taxonomy(x)))

data.to_csv(path_or_buf = output, sep = '\t', index = False)   

print('Everything done. The final output is a file named \'' + str(output) + '\' with five columns.')
print('1: contig names')
print('2: most frequent taxonomy in the contig')
print('3: taxonomy common to all proteins on the contig')
print('4: the proteins on this contig share the same taxonomy? If you find a \"Yes\" in this column, then column 2 and 3 will be the same.')
print('5: the taxonomy level of the most frequent taxonomy on the contig. Es: if the most frequent taxonomy is "Penicillium", in this column you will find "genus"')

os.remove(tempfile)


#########################################
#Script 2
"""
Step 2:
python prepare_edgeR_input_taxonomy_v2.py -o edgeR_input_taxonomy_1.contig.order.txt -c NTC_G210.count_tab.matrix -t get_contig_splitted_2_1.txt -l order
(-l you can choose a specific rank which you're interested in )
"""

count_table = BaseFolder + input2

taxonomy_file = BaseFolder + "get_contig_taxonomy_ranked.txt"
output = BaseFolder + "edgeR_input_taxonomy_1.contig.order.txt" 
level = "order"

# Prepare dataframes and dictionary
#Prepare taxonomy dictionary
ncbi = NCBITaxa()

data = pd.read_csv(taxonomy_file, sep='\t')
data.set_index(data[data.columns[0]], inplace = True)


#contigs_list = data['Contig'].values.tolist()
data["get_level"] = data.apply(lambda x: get_level_taxonomy(x["most_frequent_taxonomy"], level), axis = 1)

dic = {i[0] : i[1] for i in data[["Contig", "get_level"]].values.tolist()}

# Prepare count table
ctable = pd.read_csv(count_table, sep = '\t')
ctable.set_index(ctable['ID'], inplace = True)
ctable.drop(['ID'], axis = 1, inplace = True)
try:
    ctable.drop(['Len'], axis = 1, inplace = True)
except KeyError:
    print('')
 
    
print("Proceed to analysis")
    
#prepare empy output dataframe
outdf = pd.DataFrame(index = list(set(list(dic.values()))), columns = ctable.columns, data = 0.0)

# Fill output dataframe
for cont in dic.keys():
    for col in outdf.columns:
        outdf[col][dic[cont]] += ctable[col][cont]
  
outdf.to_csv(path_or_buf = output, sep = '\t')

#############################
##Script3

dictionary = inputUltimate_fixed#perl output
output = "edgeR_input_1.contig.txt"

# Prepare dataframes and dictionary
#Prepare dictionary from protein to EC code

dic = {}
#with open(dictionary) as t:
for line in dictionary:
    line = line.rstrip().split('\t')
    for el in line[1].rstrip().split(' '):
        if line[0] in dic:
            dic[line[0]].append(el)
        else:
            dic[line[0]] = [el]
for id in list(dic.keys()):
    dic[id] = list(set(dic[id]))


# Prepare count table
ctable = pd.read_csv(count_table, sep = '\t')
ctable.set_index(ctable['ID'], inplace = True)
ctable.drop(['ID'], axis = 1, inplace = True)
try:
    ctable.drop(['Len'], axis = 1, inplace = True)
except KeyError:
    print('')
 
print("Proceed to genewise analysis")
newindex = []
for el in dic.values():
    for code in el:
        newindex.append(code)
       
newindex = list(set(newindex))
outdf = pd.DataFrame(index = newindex, columns = ctable.columns, data = 0.0)
for gene in dic.keys():
    for col in outdf.columns:
        try:
            for code in dic[gene]:
                outdf[col][code] += ctable[col][gene]
        except KeyError:
            continue
       
outdf.to_csv(path_or_buf = BaseFolder + output, sep = '\t')

#############################
##Script4
# Define file names
#count_table = args.count_table
#dictionary = args.dictionary
level = "genus"
#taxonomy_file = args.taxonomy_file
output= BaseFolder+"ECcodes_genus_Countable.contig.txt"
proteins = BaseFolder+input1 
#tempfile = args.temp

# Prepare dataframes and dictionary
#Prepare dictionary from protein to EC code

dic = {}
for line in dictionary:
    line = line.rstrip().split('\t')
    for el in line[1].rstrip().split(' '):
        if el in dic:
            dic[el].append(line[0])
        else:
            dic[el] = [line[0]]
for id in list(dic.keys()):
    dic[id] = list(set(dic[id]))

# Prepare count table

try:
    ctable.drop(['Len'], axis = 1, inplace = True)
except KeyError:
    print('')

print("Prepare temporary files")      

## prepare protein data
with open(proteins) as f:
    with open(tempfile, 'w') as t:
        data = f.read()
        data = re.sub(r'([0-9][0-9]*_[0-9][0-9]*)_([0-9][0-9]*\t)', r'\1\t\2', data)
        data = re.sub(' ', '_', data)
        t.write(data)

annotation = pd.read_csv(tempfile, names=['contig', 'protein', 'annotation'], sep = '\t')
annotation['protein'] = annotation['protein'].astype(str)
annotation['protein'] = annotation['contig'] + '_' + annotation['protein']
   


#Prepare taxonomy of desired level
ncbi = NCBITaxa()

#read taxonomy file and prepare taxonomy dictionary of desired level
taxa = pd.read_csv(taxonomy_file, sep='\t')
taxa.set_index(taxa['Contig'], inplace = True)

dictax = {}
for cont in taxa['Contig'].tolist():
    dictax[cont] = get_level_taxonomy(taxa['most_frequent_taxonomy'][cont], level)
    
print("Proceed to contigwise analysis")

codes = []
for code in dic.keys():
        codes.append(code)
           
outdf = pd.DataFrame(index = [], columns = ctable.columns, data = 0.0)
samples = ctable.columns
outdf['EC'] = ''
outdf['taxa'] = ''

# Fill output dataframe
for code in codes:
    templist = dic[code]
    tempdf = annotation.where(annotation['protein'].isin(templist)).dropna()
    contlist = tempdf['contig'].drop_duplicates().tolist()
    for cont in contlist:
        try:
            contdf = tempdf.where(tempdf['contig'] == cont).dropna()
            protnum = len(contdf)
            tax = dictax[cont]
            subdf = outdf.where((outdf['EC'] == code) & (outdf['taxa'] == tax)).dropna()
            if len(subdf) > 1:
                print('Something unexpected happened while filling the database. Exiting now.')
                sys.exit()
            elif len(subdf) == 1:
                ind = subdf.index[0]
                for sample in samples:
                    outdf[sample][ind] += round((protnum * ctable[sample][cont]), 2)
            elif len(subdf) == 0 :
                ind = code + '___' + tax
                outdf = outdf.append(pd.DataFrame(index = [ind], columns = outdf.columns, data = 0.0))
                for sample in samples:
                    outdf[sample][ind] += round((protnum * ctable[sample][cont]), 2)
                outdf['EC'][ind] = code
                outdf['taxa'][ind] = tax
        except:
           print("error in parsing raw")

os.remove(tempfile)

cols = outdf.columns.tolist()
cols = cols[-2:] + cols[:-2]
outdf = outdf[cols]

outdf.to_csv(path_or_buf = output, sep = '\t', index = False)


#R script
data = pd.read_csv(output, sep = '\t')
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





















