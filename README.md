# VEP_Results_Parser_Merger
Parses VEP txt output (simple analysis)
by M. Joseph Tomlinson IV

Note: Code was Developed in the Abasht Laboratory at the University of 
Delaware under the supervision of Dr. Behnam Abasht website: "http://canr.udel.edu/faculty/behnam-abasht/"

Code takes in a txt file produced from an Ensembl VEP analysis
https://useast.ensembl.org/Tools/VEP

Latest version of program can handle two different types of txt files outputted by VEP without user involvement (code is smarter now)

This is a simple code that just pulls out specific columns of information and tallies the results. There are two parts of the code,
part I tallies genes and part II examines specific features of the consequence column.

#Part I. (Tallying Genes)
The code tallies the SYMBOL, Gene and Feature columns of the data, so the number of times of a 'gene' etc. shows up in the results 
is identified and allows for future pathway enrichment analysis of the final lists. 

#Part II. (Examining Consequences)
The rest of the program examines the consequences of the data (totals all features) and also examines the overall length of the 
upstream and downstream features of the data. This analysis is extremely important because the data we were examining was produced
from mRNA and so the identification of a large number of downstream_gene_variants was slightly suprising at first. However, 3'UTR can
extend very far from a gene and so the current annotation may not take this into account and so identfying the overall behavior
of the 3'UTR was extremely important. The results for the upstream and downstream are binned.

#Part III. (Merge Results)
Takes the VEP output and merges specific annotation information from the variant back into the original VCF file analyzed, so a annotated 
vcf file is now produced. Note: due to the type of merging that occurs (adding extra columns), the VCF file NO longer follows standard VCF
file format. This was specifically done to retain the original information in the VCF file. 

Required Files for Running
1. VEP txt output file
2. Original VCF file submitted to VEP
2. VRPM_Parameter_File.txt (change input file paramters)

#Running Program:
First change the parameter file's parameters. Very simple just update to your file name on the specific parameter line for the corresponding files.

Make sure to submit the correct file based on its type in the parameter file. 
(Example    --File_XYZ_Input your_file) ---note space seperated from type of input

Open up the program and run in python or double click on program to automatically run. IMPORTANT all files must be in the same directory as the program and output will occur in the same directory too. 

#Output Files
1. genes.txt (tally of genes- Symbol column)
2. protein_IDs.txt (tally of Ensemble protein IDs column)
3. ensemble_gene_IDs.txt (tally of ensemble gene IDs column)
4. upstream_variants_LOG.txt (all upstream variants in input file)
5. downstream_variants_LOG.txt (all downstream variants in input file)
6. summary_report.txt (tallying file of upstream and downstream variants with binning of distances from gene)
7. annotated_your_file_name.txt (annotated vcf file originally submitted to VEP)

