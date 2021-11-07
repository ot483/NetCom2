# 
import argparse as ap
import pandas as pd
import sys
from ete3 import NCBITaxa



#arguments
parser = ap.ArgumentParser(description='Obtain count table in which each element of desired taxonomy level will have a count equal to the sum of the counts of the contigs associated to that taxonomy.')
parser.add_argument("-t", "--taxonomy_file", help="contig taxonomy file.")
parser.add_argument("-o", "--output", help="output file")
parser.add_argument("-c", "--count_table", help="contig count table file")
parser.add_argument("-l", "--level", help="desired taxonomy level")
args = parser.parse_args()

#Check presence of arguments

if not args.count_table:
    print('A count table file is needed. Please specify it with the -i option.')
    sys.exit()
    
if not args.taxonomy_file:
    print('A dictionary file correlating contigs and taxonomy of desired level is needed. Please specify it with the -i option.')
    sys.exit()

if not args.output:
    print('A name for the output file is needed. Please specify it with the -o option.')
    sys.exit()
    
if not args.level:
    print('You need to specify the required level of taxonomy (genus, order, ecc...) with the -l option.')
    sys.exit()

# Define file names

count_table = args.count_table
taxonomy_file = args.taxonomy_file
output= args.output
level = args.level



# Prepare dataframes and dictionary

#Prepare taxonomy dictionary
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
  

ncbi = NCBITaxa()

data = pd.read_csv(taxonomy_file, sep='\t')
data.set_index(data[data.columns[0]], inplace = True)

dic = {}
for cont in data[data.columns[0]].tolist():
    dic[cont] = get_level_taxonomy(data[data.columns[1]][cont], level)
        


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

       

outdf.to_csv(path_or_buf = './' + output, sep = '\t')