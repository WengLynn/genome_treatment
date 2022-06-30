# genome_treatment: Virus Variant Identifier(VVI)

This is a R pacakge including some tools to treat aligned genome sequence.

To acess them, you can try 
   
     library(devtools)
     install_github("wuaipinglab/genome_treatment/findout")
     library(findout)
     ?findout::findout_mutation()

# 1. findout 

This is developed for SARS-CoV-2.
 
Now the available package is: 
 
#   1.1 findout_mutation (depending on tidyr, stringr)
      
   Description: When you try to call INDEL or SNP from a very large and already aligned geonome sequences on your local system, I will recommend you to use it.
      To use it try: 
       
      findout_mutation(filename, key_pattern)
     
# Parameters:

filename: The file name of your aligned sequences.

key_pattern: An unique pattern in your refference sequence in the above aligned file, such as a EPI id.

As a quality control, sequences will be removed with more than 15 "N" or 50 other degenerate bases like "R".

#    1.2 NTtoAA (depending on ape, dplyr)

     Description: If you have la list of nucleotide mutation and want to figure out which protein located on or whether lead to a non-synonymous mutation or not, I would like to introduce this package to you.
     
     To use it try: 
       
      NTtoAA(Seq_fasta, ORF_table_file, SNP_table_file)
      
# Parameters:

Seq_fasta: reference sequence, make sure this seqeunce has the closest evolution distance to your input.

ORF_table_file: Address to a csv file list a Open Reading Frame Table to your reference sequence.

SNP_table_file: Address to a csv file list the SNP intresting you, such as A21222T. It should be in nucleotide.

ORFStartCol = "Start": The colname of the Start site in ORF_table_file

ORFEndCol = "End": The colname of the End site in ORF_table_file

ORFGeneCol="Gene": The colname of the Gene name in ORF_table_file

SNP = "SNP": The colname of the SNP in SNP_table_file

# Good Luck


