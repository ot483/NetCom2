# Microbiome_2021
Supplementary materials for Berihu et al. \
The input/output files used/generated with code below are all available on Drive: https://volcanicenter-my.sharepoint.com/:f:/g/personal/ofirt_volcani_agri_gov_il/EsPau_uqpolHk37VuXWMqMIB0J-Ey1I-Pstl2gt0k48G8A
# Code for merging MEGAN annotations tables with count tables
A count table was created using BWA mapping software. Files for  both functional and taxonomic annotations of genes were independently generated with MEGAN software. The count data are presented as a table which reports, for each sample, the number of sequence fragments that have been assigned to each contigs or genes. To merge count table with functional and taxonomic keys, we wrote the following codes written in python: \
Step 1: "get_contig_taxonomy_v3.py" takes as input megan taxonomic annotation for each gene and defines the frequent taxonomy of genes associated with a contig and prints the taxonomy rank of the contig. \
Step 2: "Prepare_edgeR_input_taxonomy_v2.py" merges the count table and the taxonomic annotations of the contigs (generated in step 1), creating a count table based on the selected taxonomic level. \
Step 3: "script_1kegg_parser.pl" filters only contigs/genes with EC/KEGG identifiers from MEGAN functional file: "ReadName_to_KeggName". \
Step 4: "prepare_edgeR_input_v2.py" merges the count table, MEGAN file with functional annotation, and the output file from Step 3 creating a count table based on the functional keys (EC/KEGG). \
Step 5: "prepare_ECcodes_taxonomy_ctable_contigwise_v2.py" merges a count table with both taxonomic and functional key for selected scheme/taxonomic levels. \
Step 6: "knockout-ecs_interactive_modified.R" uses the output file generated in Step 5 and generates a diversity score for the taxonomic distribution of each enzyme. \
All input and output files are available on the drive https://volcanicenter-my.sharepoint.com/:f:/g/personal/ofirt_volcani_agri_gov_il/EsPau_uqpolHk37VuXWMqMIB0J-Ey1I-Pstl2gt0k48G8A
# Using edgeR to obtain a list of significantly different contigs/genes using_edgeR_generic_2.R
The edge R scripts identify differentially abundant taxonomic/functional groups and requires six parameters: \
1. Count table with a selected key index. 
2. Metadata file, with at least two columns: one for the samples ID and one assigning each sample to a group. 
3. Name of output files (the program generates multiple outputs,  all with unique prefix). 
4. Name of the metadata file column dividing the samples into the groups of replicates. 
5. Text file with the treatments being compared (e.g., NTC G210 vs NTC M26). The names of the selected treatment should be separated by ' = ' delimiter. \
6. FDR threshold used to determine significance. For our analysis, we used 0.05 . 
Example for running the code: \
Rscript using_edgeR_generic.R. Count_Table_Order_NTC.txt sample_metadata.txt Order_0.05_ Tretment_Rootstock contrasts_used.txt 0.05. 
The script will produce, for each treatment, a table file with the results of the differential analysis and a plot showing the distribution of differentially represented elements. \
Network analysis Code for network construction and visualization was deposited in https://github.com/ot483/NetCom
# Removal Network
The script "Removal_Network_order_1_1.24.py" was used as a reference for community 'knockouts' simulations in which selected taxonomic groups were removed. In each of the removal iterations, all edges -enzymes representing metabolic functions, specifically dominated by a taxonomic group, taxa-dominated enzymes were removed from the original enzyme set. The impact of the removal group was estimated according to differences in the number of metabolites between the network expanded from the truncated enzyme set, and the reference meta-network.
The four required DB files are available in the Scripts folder.
