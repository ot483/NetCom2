# 
import argparse as ap
import pandas as pd
import sys
from ete3 import NCBITaxa
import os
import re



#arguments

parser = ap.ArgumentParser(description='Obtain count table in which each row is a combination of a function and a taxonomy level (genus, orders, ecc...')
parser.add_argument("-r", "--temp", help="temporary file. Only necessary if -g option is not used.")
parser.add_argument("-d", "--dictionary", help="dictionary file correlating proteins and functions (EC codes or KEGG names)")
parser.add_argument("-t", "--taxonomy_file", help="contig taxonomy file.")
parser.add_argument("-o", "--output", help="output file")
parser.add_argument("-p", "--proteins", help="file correaling each protein to its contig, like the DIAMOND output file.")
parser.add_argument("-l", "--level", help="desired taxonomy level")
parser.add_argument("-c", "--count_table", help="count table file")
args = parser.parse_args()

#Check presence of arguments


if not args.dictionary:
    print('A dictionary file correlating proteins and functions is needed. Please specify it with the -d option.')
    sys.exit()

if not args.level:
    print('A target taxonomy level (es: genus) is needed. Please specify it with the -l option.')
    sys.exit()

if not args.taxonomy_file:
    print('A contig taxonomy file is needed. Please specify it with the -t option.')
    sys.exit()

if not args.output:
    print('A name for the output file is needed. Please specify it with the -o option.')
    sys.exit()

if not args.proteins:
    print('A protein file is needed. Please specify it with the -i option.')
    sys.exit()
    
if not args.temp:
    print('A name for the temporary file is needed. Please specify it with the -r option.')
    sys.exit()
    
# Define file names


count_table = args.count_table
dictionary = args.dictionary
level = args.level
taxonomy_file = args.taxonomy_file
output= args.output
proteins = args.proteins
tempfile = args.temp




# Prepare dataframes and dictionary
#Prepare dictionary from protein to EC code

dic = {}
with open(dictionary) as t:
    for line in t:
        line = line.rstrip().split('\t')
        for el in line[1].rstrip().split(' '):
            if el in dic:
                dic[el].append(line[0])
            else:
                dic[el] = [line[0]]
    for id in list(dic.keys()):
        dic[id] = list(set(dic[id]))

# Prepare count table

ctable = pd.read_csv(count_table, sep = '\t')
ctable.set_index(ctable[ctable.columns[0]], inplace = True)
ctable.drop([ctable.columns[0]], axis = 1, inplace = True)
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
    #Define functions
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

os.remove(tempfile)

cols = outdf.columns.tolist()
cols = cols[-2:] + cols[:-2]
outdf = outdf[cols]


outdf.to_csv(path_or_buf = './' + output, sep = '\t', index = False)

