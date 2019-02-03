# Divergence-based Introgression Polarization (DIP)

![](images/triple-DIP.png)

Introgressive hybridization is an impactful evolutionary process. 
While there are many tools for identifying the presence of introgression, few explicitly determine the direction in which genetic material was transferred during introgression.
DIP is a tool for polarizing both unidirectional and asymmetric bidirectional introgression to determine the donor and recipient taxa.
DIP works in a widely-used four-taxon context and is built to analyze either whole-chromosome alignments or single-locus alignments spanning the genome.

## DIP References:
*_<DIP paper, link coming soon>_*

*Forsythe ES, Nelson AD, Beilstein MA. Biased gene retention in the face of massive nuclear introgression obscures species relationships (Submitted). bioRxiv. Available from: https://www.biorxiv.org/content/early/2018/10/18/197087?%3Fcollection=*


## Contents of this repository
This repository contains scripts for performing Divergence-based Introgression Polarization (DIP) analyses.

### *DIP.R*

DIP.R is used to perform 1x, 2x, and 3x-DIP analyses. 
DIP.R is called from the command line using the rscript command.
The input data are single-gene multiple sequence alignments in fasta or phylip format.
DIP was designed for and tested on datasets composed of roughly 1000-10,000 single locus alignments.
Alignment files should be stored be stored in the same directory. 
The path to this directory specified by the user when DIP.R is called.
The names of all alignment file should contain a common string, that is also specified by the user when DIP.R is called.
The taxon identifiers/sequence names should contain unique strings specifying the species name (e.g. SpeciesX_000000001, SpeciesY_000000001, etc...)

*Running DIP.R*

DIP.R is called from the command line with nine arguments as follows:

`Rscript --vanilla DIP.R <A1:Jobname> <A2:type_of_alignment> <A3:directory_containing_alignments_and_listfile> <A4:Simlulation replicate> <A5:Total_number_of_taxa_in_alignments> <A6:P1_taxon_string> <A7:P2_taxon_string> <A8:P3_taxon_string> <A9:Outgroup_taxon_string>`

+ *<A1:Jobname>*

This should be used to specify a string that is common to all alignment file names

+ *A2:type_of_alignment*

The alignment type. Either "fasta" or "phylip"

+ *A3:directory_containing_alignments_and_listfile*

The full path to the directory containing the single-locus alignments. End with a "/"

+ *<A4:Simulation replicate>*

This is only used for simulation analyses where the genome alignment step was replicated. For analysis of empirical data, make this the exact same string used in A1

+ *A5:Total_number_of_taxa_in_alignments*

The total number of taxa in each alignment. If this number varies, put the minimum number of taxa expected. Should be at least four.

+ *A6:P1_taxon_string*

A string found in the sequence identifiers of the species used as P1

+ *A7:P2_taxon_string*

A string found in the sequence identifiers of the species used as P2

+ *A8:P3_taxon_string*

A string found in the sequence identifiers of the species used as P3

+ *A9:Outgroup_taxon_string*

A string found in the sequence identifiers of the species used as the outgroup









