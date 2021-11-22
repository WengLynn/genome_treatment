# genome_treatment
This is a R pacakge including some tools to treat aligned genome sequence.
To acess them, you can try 
  install_github("WengLynn/genome_treatment/package_name")
Now the available package is: 
  1. findout (This is developed for SARS-CoV-2.)
      1.1 findout_mutation (depending on tidyr, stringr)
      Description: When you try to call INDEL or SNP from a very large and already aligned geonome sequences on your local system, I will recommend you to use it.
      To use it try: findout_mutation(filename,key_pattern)
      
      1.2 findout_whichprotein (coming soon)
      Description: When you try to know which protein or protreins your INDEL or SNP located in from a large table, I will recommend you to use it.
      
      1.3 findout_overlap (coming soon)
      Description: When you try to stick your INDEL region to get their overlap region from a large table, I will recommend you to use it.
      
      1.4 findout_whichsite (coming soon)
      Description: When you try to now where your sites aganist your reference sequence, I will recommend you to use it.
      
      1.5 findout_AAchange (coming soon)
      Description: When you try to the impact of your INDEL or SNP on amino acid change, I will recommend you to use it.
      
      1.6 findout_mapacross (coming soon)
      Description: When you try to map your INDEL or SNP from SARS to SARS2, I will recommend you to use it.
      
      ..ongoing..
       
  2.  shorten (This is developed NOT for SARS-CoV-2.)
      2.1 shorten_site
      Description: When you try to remove your conserved sites to reduce sequence length for further analysis such as calling adaptive mutations, I will recommend you to use it.
  
     ..ongoing..
     
Slow Down and Good Night
