<img src="dep_sign.png" width=120, height=120 align="left" />

# NetCom2

A pipeline for reproducing sequence processing and analysis of metagenomics data produced for Berihu et al ("A framework for the targeted recruitment of crop-beneficial soil taxa based on network analysis of metagenomics data")

## Benchmarks

All datasets mentioned in text are available in **https://volcanicenter-my.sharepoint.com/:f:/g/personal/shmedina_volcani_agri_gov_il/EicqyHpmSmBJmwgHO0oqJ8MBcxNtvXw0KUs6hYyyQ4wvzA?e=khRLmH**.

## Dependencies

* [python (version >= 3.8)]
* [matplotlib (version 3.3.3)]
* [networkx (version 2.5)]
* [pandas (version 1.1.4)]
* [skbio (version 0.5.6)]
* [numpy (version 1.19.4)]
* [ete3 (version 3.1.2)]

## Installation
### Download and install NetCom2 on linux (Ubuntu 20.04)
### Create virtual environment and install dependencies

```shell
# Create venv #
mkdir netcom2
cd netcom2
virtualenv netcom2
source netcom2/bin/activate
```

### Download and install NetCom2 on Windows
```shell
mkdir netcom2
cd netcom2
python -m venv C:\path\netcom2
C:\path\netcom2\Scripts\Activate
```

```shell
# You may need to install svn and pip3 #
sudo apt install subversion
svn export https://github.com/ot483/NetCom2/trunk/Scripts
curl -L https://raw.githubusercontent.com/ot483/NetCom2/main/requirements.txt
pip install -r requirements.txt
```



## Tutorial

### Run NetCom2

#### Get contig taxonomy

The script determines the taxonomic classification of each contig according to the annotations of its respective genes. The script takes as input a Megan output file with taxonomic annotation for each gene and determines the taxonomic annotation of the corresponding contigs according to most frequent taxonomy of the associated genes.

Input: **NTC_G210-ReadName_to_TaxonPath.txt**

Output: **get_contig_taxonomy_ranked_NTC_G210.txt**

```shell
python get_contig_taxonomy.py -i NTC_G210-ReadName_to_TaxonPath.txt -t temp_1 -o get_contig_taxonomy_ranked_NTC_G210.txt --get_rank
```

#### Prepare taxonomy count table

The script merges the count table (constructed for contigs, Methods) and the taxonomic annotations of the contigs (generated in step 1) and creates a count table based on the selected taxonomic level by merging respective contigs.

Input: **get_contig_taxonomy_ranked_NTC_G210.txt**, **NTC_G210.count_tab.matrix**

Output: **count_table_taxonomy_NTC_G210.contig.genus.txt**

```shell
python prepare _taxonomy_count_table.py -o count_table_taxonomy_NTC_G210.contig.genus.txt -c NTC_G210.count_tab.matrix -t get_contig_taxonomy_ranked_NTC_G210.txt -l genus
```
-l : specifying taxonomic rank according to ncbi taxonomy etc species,genus,order,phylum

#### Filter enzymatic functions

The script takes as input a Megan output file with KEGG functional annotations and filters only entities with EC/KEGGs KO identifiers.

Input: **NTC_G210-ReadName_to_KeggName_Ultimate.txt**

Output: **protein_EC_codes_total_NTC_G210.txt**

```shell
Perl Filter_enzymatic_functions.pl
Change parameters: my $file ="NTC_G210-ReadName_to_KeggName_Ultimate"
open IN ,$file
open (OUT,">protein_EC_codes_total_NTC_G210.txt")
```

#### Prepare function count table

The script merges the count table with functional annotations retrieved from MEGAN and filtered in Step 3, creating a count table based on the functional keys (EC/KEGGs KO).

Input: **protein_EC_codes_total_NTC_G210.txt**, **NTC_G210.count_tab.matrix**, **NTC_G210-ReadName_to_TaxonPath.txt**

Output: **count_table_function_NTC_G210.contig.txt**

```shell
python prepare_function_count_table.py -d protein_EC_codes_total_NTC_G210.txt -t temp_1 -o count_table_function NTC_G210.contig.txt -c NTC_G210.count_tab.matrix -p NTC_G210-ReadName_to_TaxonPath.txt
```

#### Prepare function taxonomy count table contigwise

The script merges the count table according to a taxonomic (Step 2) and functional (Step 4) keys into a double-key (taxonomic/functional) count table. The output file details the taxonomic distribution of reads assigned to each enzymes at a selected taxonomic level.

Input: **protein_EC_codes_total_NTC_G210.txt**, **NTC_G210.count_tab.matrix**, **NTC_G210-ReadName_to_TaxonPath.txt**, **get_contig_taxonomy_NTC_G210.txt**

Output: **function_taxonomy_countable_NTC_G210.contig.genus.txt**

```shell
python prepare_function_taxonomy_ctable_contigwise.py -d protein_EC_codes_total_NTC_G210.txt -r temp_1 –o function_taxonomy _countable_NTC_G210.contig.genus.txt -c NTC_G210.count_tab.matrix -p NTC_G210-ReadName_to_TaxonPath.txt -t get_contig_taxonomy_NTC_G210.txt -l genus
```
-l : specifying taxonomic rank according to ncbi taxonomy etc species,genus,order,phylum


#### Calculate enzyme diversity score

The script takes as input the enzyme – diversity table generated at Step 5 and calculates a diversity score for each enzyme, indicating whether it is dominantly present in a specific taxonomic group (at a selected taxonomic level) or widely distributed across microbial community.

Input: **function_taxonomy_countable_BjSa_NTC.contig.genus.txt**

Output: **ECs_with_dominant_taxon.txt**, **ECs_with_dominant_taxon_knockout.txt**, **RscriptOutput.csv**, **RscriptOutput_filter_IdentMostFreq.csv**, **RscriptOutput_filter_IdentMostFreqbyEnzyme.csv**, **RscriptOutput_filter_dominance.csv**

```shell
python calculate_enzyme_diveristy_score.py function_taxonomy_countable_BjSa_NTC.contig.genus.txt /path/to/input/file/
```

### Above pipeline implemented as a single script using multi-thread processing.

A pipeline carrying out together steps 1-6 above

Input: **NTC_G210-ReadName_to_TaxonPath.txt**, **NTC_G210.count_tab.matrix**, **NTC_G210-ReadName_to_KeggName_Ultimate.txt**

Output: output of steps 1-6

Three input files are required - ReadName_to_TaxonPath.txt, count_tab.matrix and ReadName_to_KeggName_Ultimate.txt"
There are 4 arguments that have to be stated by the following order: ReadName_to_TaxonPath.txt
count_tab.matrix
ReadName_to_KeggName_Ultimate
BaseFolder - where the script and input files are. python 

```shell
Removal_Pipeline.py ReadName_to_TaxonPath.txt count_tab.matrix ReadName_to_KeggName_Ultimate /Path/To/InputFiles/
```

#### Determine differential abundance

The script determines differential abundance between entities whose relative abundance is described in a count table. The script will produce, for each treatment, a table file with the results of the differential analysis and a plot showing the distribution of differentially represented elements.

Input: **count_table_taxonomy_genus_BjSaVSNTC.txt**, **sample_metadata.txt** (Metadata file, with at least four columns: samples ID, name of the rootstock, name of treatment and name of rootstock treatment), **contrasts_used_BjSa_NTC.txt**

Output: analysis results

```shell
Rscript using_edgeR_generic.R. count_table_taxonomy_genus_ BjSaVSNTC.txt sample_metadata.txt Perfix_for_output_files Treatment_Rootstock contrasts_used_BjSa_NTC.txt 0.05
```
Parameters: count table, metadata file, prefix for output files, compared treatments as specified by column names in metadata file, text file with the treatments being compared (e.g., NTC vs BjSa), FDR threshold.
0.05-FDR threshold used to determine significance.


The script tool for predicting metabolic activities in microbial communities based on network-based interpretation of assembled and annotated metagenomics data. The algorithm takes as input an EdgeR output file that provides information on the differential abundance of enzymatic reactions in two different treatments (Step 8). Enzymes are classified as associated with Treatment_1, Treatment_2 or not associated. The algorithm generates as output: (i) Lists of differentially abundant enzymes and their pathway association. (ii) Prediction of environmental resources that are unique to each treatment and their pathway association. (iii) Prediction of environmental compounds that are produced by the microbial community and pathway association of compounds that are treatment-specific. (iv) Network visualization of enzymes, environmental resources and produced compounds that are treatment specific (2 & 3D).

Code and instructions are available at https://github.com/ot483/NetCom

#### Community 'knockouts'  simulations

The script conducts community 'knockouts' simulations in which selected taxonomic groups are removed from the community network by eliminating enzymes associated with the group (according to the scores in step 6). The impact of the removal group is estimated according to differences in the number of metabolites between the network expanded from the truncated enzyme set, and the reference meta-network that includes the full set of enzymatic functions

Input: Knock_out_file_step_6_update_name, **Env.txt**, dictionary files (From the NetCom package)

Output: analysis results

```shell
svn export https://github.com/ot483/NetCom2/trunk/Dict
Python Sim4RemovalNet.py
```
Input/output files names can be updated in the script.
Knock_out_file_step_6_update_name is the output of step 6.
Env.txt is an output of step 9 – product (ii) Prediction of environmental resources that are unique to each treatment and their pathway association. The script requires dictionary files provided in Dict folder: **compounds_lables_jun_1.txt**, **ec_reac_mapping_jun.txt**, **full_enzymes_labels_jun.txt**, **reactions_3_balanced.txt**

#### Visualization

The script takes as input lists of enzymes, predicted environmental resources and unique compounds and produces a network representation.

Input: **allCompounds_BjSa.txt**, **allCompounds_NTC.txt**, **Compounds_BjSa_Order.txt**, **EC_ALL.txt, ECs_BjSa.txt**, **ECs_NTC.txt, enzymes_BjSa_order.txt**, **pathways_BjSa_order.txt**, **seeds_BjSa.txt**, **seeds_NTC.txt**, **compounds_lables_jun_1.txt**, **ec_reac_mapping_jun.txt**, **full_enzymes_labels_jun.txt**, **reactions_3_balanced.txt**

Output: **BjSa_order_removal_network.pdf**, **BjSa_oder_removal_network.png**

The user should change the file names in the script - lines 80-128 

```shell
python3 Network_Figures.py
```

OR

The script can be used by arguments as below, just uncomment lines 33-78 and comment lines 80-128


```shell
python3 Network_Figures.py <This is base folder> <drop fragments with size smaller then> \
                          <integer of hubness to filter> <prefix> <All_compounds_B_input txt file> \
                          <Seeds_B_input txt file> <ECs_All_input txt file> \
                          <patches_Compounds> <patches_Enzymes> <Pathways_Compounds>
```



## References

