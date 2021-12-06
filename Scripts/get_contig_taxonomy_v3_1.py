# 
import re
import argparse as ap
import pandas as pd
import sys
import os
import numpy as np

#arguments
parser = ap.ArgumentParser(description='Get taxonomy of contigs from DIAMOND output')
parser.add_argument("-i", "--input", help="input file")
parser.add_argument("-t", "--temp", help="temporary file")
parser.add_argument("-o", "--output", help="output file")
parser.add_argument("--get_rank", help="get from ncbi the rank of contig taxonomy lavels. It greatly increases run time, and it requires the ete3 module.", action="store_true")
args = parser.parse_args()

#Check presence of arguments

if not args.input:
    print('An input file is needed. Please specify it with the -i option.')
    sys.exit()

if not args.temp:
    print('A name for the temporary file is needed. Please specify it with the -t option.')
    sys.exit()

if not args.output:
    print('A name for the output file is needed. Please specify it with the -o option.')
    sys.exit()

# Define file names

file = args.input
tempfile = args.temp
output= args.output

 



#Check input file
try:
   with open(file) as f:
       print('I will now start to work. Feel free to do something relaxing while you wait.')
except FileNotFoundError:
    print('It seems that the input file you specified does not exist. I am very sorry, but I cannot proceed :(.')
    sys.exit()

#At first I create the temporary file with the data formatted as I need
def create_temp_file(input_, tempfile):
    with open(input_) as f:
        with open(tempfile, 'w') as t:
            data = f.read()
            data = re.sub(r'([0-9][0-9]*_[0-9][0-9]*)_([0-9][0-9]*\t)', r'\1\t\2', data)
            data = re.sub(' ', '_', data)
            t.write(data)



data = pd.read_csv(tempfile, names=['contig', 'protein', 'annotation'], sep = '\t')
contlist = set(data['contig'].tolist())
contnum = len(contlist)

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

        
            
create_temp_file(file, tempfile)


from multiprocessing import Pool

#number_current_contig = 0

#o = open(output, 'w')
#o.write('Contig' + '\t' + 'most_frequent_taxonomy' + '\t' + 'common_taxonomy' + '\t' + 'Same_taxonomy?' + '\n')
#results_list.append(['Contig' + '\t' + 'most_frequent_taxonomy' + '\t' + 'common_taxonomy' + '\t' + 'Same_taxonomy?' + '\n'])

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


results_list_final = Pool().map(manipulate, list(contlist))

results_list_final_flatten = [ent for sublist in results_list_final for ent in sublist]
results_list_final_flatten = [['Contig' + '\t' + 'most_frequent_taxonomy' + '\t' + 'common_taxonomy' + '\t' + 'Same_taxonomy?' + '\n']]+ results_list_final_flatten

o = open(output, "w")
for element in results_list_final_flatten:    
    o.write(element[0])
o.close()


from ete3 import NCBITaxa

ncbi = NCBITaxa()
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