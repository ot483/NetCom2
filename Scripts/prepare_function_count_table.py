# 
import re
import argparse as ap
import pandas as pd
import sys
import os




#arguments
parser = ap.ArgumentParser(description='Obtain df in which each EC code has, for each protein referring to it, a count value equal to the reads mapped to the contig with the protein. This will result in a bias towards large contigs, unless the count table is already normalized according to contig size.')
parser.add_argument("-d", "--dictionary", help="dictionary file correlating proteins and functions (EC codes or KEGG names)")
parser.add_argument("-t", "--temp", help="temporary file. Only necessary if -g option is not used.")
parser.add_argument("-o", "--output", help="output file")
parser.add_argument("-c", "--count_table", help="count table file. By default the program assumes it is based on contigs. if it is based on genes then activate the -g option")
parser.add_argument("-g", "--gene_count_table", help="this option means that the count table is based on genes and not on contigs. If used, each EC code will have a count equal to the sum of the counts of the proteins presenting the EC code.", action="store_true")
parser.add_argument("-p", "--proteins", help="file correaling each protein to its contig, like the DIAMOND output file. Only necessary if -g option is not used.")
args = parser.parse_args()

#Check presence of arguments

if not args.count_table:
    print('A count table file is needed. Please specify it with the -i option.')
    sys.exit()
    
if not args.dictionary:
    print('A dictionary file correlating proteins and functions is needed. Please specify it with the -i option.')
    sys.exit()

if not args.output:
    print('A name for the output file is needed. Please specify it with the -o option.')
    sys.exit()

# Define file names

count_table = args.count_table
dictionary = args.dictionary
output= args.output
gene_count_table = False
if args.gene_count_table:
    gene_count_table = True


# Prepare dataframes and dictionary
#Prepare dictionary from protein to EC code

dic = {}
with open(dictionary) as t:
    for line in t:
        line = line.rstrip().split('\t')
        if gene_count_table == False:
            for el in line[1].rstrip().split(' '):
                if el in dic:
                    dic[el].append(line[0])
                else:
                    dic[el] = [line[0]]
        elif gene_count_table == True:
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
 



if gene_count_table == False:
    
    if not args.proteins:
        print('A protein file is needed. Please specify it with the -i option.')
        sys.exit()
    if not args.temp:
        print('A name for the temporary file is needed. Please specify it with the -t option.')
        sys.exit()
      
    proteins = args.proteins
    tempfile = args.temp
    
    print("Proceed to contigwise analysis")
    #prepare empy output dataframe
    outdf = pd.DataFrame(index = list(set(list(dic.keys()))), columns = ctable.columns, data = 0.0)
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

        
# Fill output dataframe
    for code in outdf.index:
        templist = dic[code]
        tempdf = annotation.where(annotation['protein'].isin(templist)).dropna()
        contlist = tempdf['contig'].drop_duplicates().tolist()
        for col in outdf.columns:
            reads = 0.0
            for cont in contlist:
                contdf = tempdf.where(tempdf['contig'] == cont).dropna()
                protnum = len(contdf)
                reads += round((protnum * ctable[col][cont]), 2)
                outdf[col][code] = reads
    
    
    os.remove(tempfile)
        
 

elif gene_count_table == True:
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
       

outdf.to_csv(path_or_buf = './' + output, sep = '\t')
