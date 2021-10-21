# BMC_2021
Supplementary materials for Berihu et al

Merging MEGAN annotations with counא table 

We created a count table using BWA mapping software, and each treatment was run on MEGAN to create functional or taxonomic keys. The count data are presented as a table which reports, for each sample, the number of sequence fragments that have been assigned to each contigs or genes. 
To merge count table with functional and taxonomic keys, we wrote the following codes on python.

The first code defines the contig/gene frequent taxonomy and prints the taxonomy rank. 
get_contig_taxonomy_v3.py


The code merges the count table and the new file created from the first code to built count table based on the taxonomic level.
prepare_edgeR_input_taxonomy_v2.py


Cut only contigs/genec with EC/KEGG fron Megan functioanal file "ReadName_to_KeggName".
script_1kegg_parser.pl


The code merges count table, MEGAN file with taxonomy annotation, and the new output file from "script_1kegg_parser.pl" to built count table based on the functional keys (EC/KEGG).
prepare_edgeR_input_v2.py


Merges count table, MEGAN file with taxonomy annotation, and the new output file from "script_1kegg_parser.pl" to built count table based on the functional keys (EC/KEGG) with taxonomic level.
prepare_ECcodes_taxonomy_ctable_contigwise_v2.py


Using edgeR to obtain a list of significantly different contigs/genes
using_edgeR_generic_2.R


The input/output files are available on Drive 

# EdgeR
The edge R script requires 6 arguments, the following order:
1: pathway to the input count table used
2: pathway to the metadata file, with at least two columns: one for the samples ID and one assigning each sampe to a group
3: base name for output files (the program generates multiple outputs, this is just the pathway to them and the "base" of the name, such as "edgeR_results"
4: name of the metadata file column dividin the samples into the groups of replicates
5: text file with the contrasts used in the analysis. It is important that the contrast names and the contrast structure (see example files) are separates by ' = '. In the example, we want to compare The conventional israelian group of samples to the organic israelian group and to the average of the american conventional groups.
6: FDR threshold used to establish significance. For our analysis, I used 0.05

Example:
Rscript using_edgeR_generic.R Count_Table_Order_NTC.txt sample_metadata.txt Order_0.05_ Tretment_Rootstock contrasts_used.txt 0.05

The script will produce, for each contrast, a table file with the results of the differential analysis and a plot showing the distribution of differentially represented elements.
