findout_mutation_NTtoAA <- function(Seq_fasta ,ORF_table_file,SNP_table_file,ORFStartCol = "Start",ORFEndCol = "End",ORFGeneCol="Gene",SNP = "SNP"){
  ORF_table = read.csv(`ORF_table_file`)
  SNP_table = read.csv(`SNP_table_file`)
  library(ape)
  library(dplyr)
  seq = read.table(Seq_fasta,sep = "\t")
  if(length(seq[,1])>2){
    Seq_fasta  = read.FASTA(Seq_fasta);write.FASTA(Seq_fasta,"temp.fasta");seq=read.table("temp.fasta",sep = "\t");system("rm temp.fasta")
  }
  seq=unlist(strsplit(seq[2,1],""))
  ORF_table[which(ORF_table[`ORFEndCol`]>ORF_table[`ORFStartCol`]),"Dir"] = "+"
  ORF_table[which(ORF_table[`ORFEndCol`]<ORF_table[`ORFStartCol`]),"Dir"] = "-"
  SNP_table = SNP_table[grep("[A-Z]",unlist(SNP_table[`SNP`])),]
  SNP_table$"loc" = gsub("[A-Z]","",unlist(SNP_table[`SNP`]))
  if(length(grep("-",SNP_table[,`SNP`]))>0){SNP_table = SNP_table[-grep("-",SNP_table[,`SNP`]),]}
  if(length(grep("_",SNP_table[,`SNP`]))>0){SNP_table = SNP_table[-grep("_",SNP_table[,`SNP`]),]}
  if(length(grep("\\*",SNP_table[,`SNP`]))>0){SNP_table = SNP_table[-grep("\\*",SNP_table[,`SNP`]),]}
  if(length(grep("\\#",SNP_table[,`SNP`]))>0){SNP_table = SNP_table[-grep("\\#",SNP_table[,`SNP`]),]}
  SNP_table$"loc" = as.numeric(SNP_table$"loc")
  SNP_table$"NNfrom" = gsub("[0-9].*","",unlist(SNP_table[`SNP`]))
  SNP_table$"NNto" = gsub(".*[0-9]","",unlist(SNP_table[`SNP`]))
  seq_to = seq
  for (i in 1:length(SNP_table[,1])) {
    seq_to[SNP_table[i,"loc"]] = SNP_table[i,"NNto"]
    if(seq_to[SNP_table[i,"loc"]] != SNP_table[i,"NNfrom"]){
      warning(paste0("Location ",SNP_table[i,"loc"],":",SNP_table[i,"NNfrom"],"in SNP_From is not match to the reference ",seq_to[SNP_table[i,"loc"]]))
    }
  }
  # anno --------------------------------------------------------------------

  for (i in 1:length(SNP_table[,1])) {
    SNP_table[which(SNP_table$loc>=as.numeric(ORF_table[i,`ORFStartCol`]) & SNP_table$loc<=as.numeric(ORF_table[i,`ORFEndCol`])),  "protein"] = ORF_table[i,`ORFGeneCol`]
    SNP_table[which(SNP_table$loc>=as.numeric(ORF_table[i,`ORFEndCol`])   & SNP_table$loc<=as.numeric(ORF_table[i,`ORFStartCol`])),"protein"] = ORF_table[i,`ORFGeneCol`]
  }
  SNP_table =left_join(SNP_table,ORF_table,by = c("protein"=`ORFGeneCol`))

  conwant = which(colnames(SNP_table)%in%c(c(`SNP`,`ORFStartCol`,`ORFEndCol`,`ORFGeneCol`,"loc","Dir","protein")))
  SNP_table = SNP_table[conwant]
  for (i in which(SNP_table$Dir == "+")){
    SNP_table[i,"shift"] = (SNP_table[i,"loc"]- SNP_table[i,`ORFStartCol`]+1)%%3
    SNP_table[i,"AAloc"] = ceiling((SNP_table[i,"loc"]- SNP_table[i,`ORFStartCol`]+1)/3)}
  for (i in which(SNP_table$Dir == "-")){
    SNP_table[i,"shift"] = (abs(SNP_table[i,"loc"]- SNP_table[i,`ORFEndCol`]-1))%%3
    SNP_table[i,"AAloc"] = ceiling((abs(SNP_table[i,"loc"]- SNP_table[i,`ORFEndCol`]-1))/3)}
  for (i in which(SNP_table$Dir == "+")){
    loc =SNP_table[i,"loc"]
    if ( SNP_table[i,"shift"] == 0){SNP_table[i,"AAfrom_codon"] = paste0(seq[loc-2],seq[loc-1],seq[loc])
    SNP_table[i,"AAto_codon"] = paste0(seq_to[loc-2],seq_to[loc-1],seq_to[loc])}
    if ( SNP_table[i,"shift"] == 2){SNP_table[i,"AAfrom_codon"] = paste0(seq[loc-1],seq[loc],seq[loc+1])
    SNP_table[i,"AAto_codon"] = paste0(seq_to[loc-1],seq_to[loc],seq_to[loc+1])}
    if ( SNP_table[i,"shift"] == 1){SNP_table[i,"AAfrom_codon"] = paste0(seq[loc],seq[loc+1],seq[loc+2])
    SNP_table[i,"AAto_codon"] = paste0(seq_to[loc],seq_to[loc+1],seq_to[loc+2])}
  }

  for (i in which(SNP_table$Dir == "-")){
    loc =SNP_table[i,"loc"]
    if ( SNP_table[i,"shift"] == 1){SNP_table[i,"AAfrom_codon"] = paste0(rev(c(seq[loc-2],seq[loc-1],seq[loc])),collapse = "")
    SNP_table[i,"AAto_codon"] = paste0(rev(c(seq_to[loc-2],seq_to[loc-1],seq_to[loc])),collapse = "")}
    if ( SNP_table[i,"shift"] == 2){SNP_table[i,"AAfrom_codon"] = paste0(rev(c(seq[loc-1],seq[loc],seq[loc+1])),collapse = "")
    SNP_table[i,"AAto_codon"] = paste0(rev(c(seq_to[loc-1],seq_to[loc],seq_to[loc+1])),collapse = "")}
    if ( SNP_table[i,"shift"] == 0){SNP_table[i,"AAfrom_codon"] = paste0(rev(c(seq[loc],seq[loc+1],seq[loc+2])),collapse = "")
    SNP_table[i,"AAto_codon"] = paste0(rev(c(seq_to[loc],seq_to[loc+1],seq_to[loc+2])),collapse = "")}
  }
  for (i in which(SNP_table$Dir == "-")){
    temp = SNP_table[i,"AAto_codon"]
    temp = unlist(strsplit(temp,""))
    temp_A = which(temp == "A")
    temp_C = which(temp == "C")
    temp_T = which(temp == "T")
    temp_G = which(temp == "G")
    temp[temp_A] = "T"
    temp[temp_T] = "A"
    temp[temp_C] = "G"
    temp[temp_G] = "C"
    temp = paste0(temp,collapse = "")
    SNP_table[i,"AAto_codon"] =temp
    temp = SNP_table[i,"AAfrom_codon"]
    temp = unlist(strsplit(temp,""))
    temp_A = which(temp == "A")
    temp_C = which(temp == "C")
    temp_T = which(temp == "T")
    temp_G = which(temp == "G")
    temp[temp_A] = "T"
    temp[temp_T] = "A"
    temp[temp_C] = "G"
    temp[temp_G] = "C"
    temp = paste0(temp,collapse = "")
    SNP_table[i,"AAfrom_codon"] =temp

  }


  AAtrans = read.delim("AAtrans_Single.txt")
  for (i in 1:length(SNP_table[,1])){
    if(!is.na(SNP_table[i,"protein"])){
      ko=c()
      kk = SNP_table[i,"AAfrom_codon"]
      for (k in 1:6) {ko=c(ko,rownames(AAtrans)[which(AAtrans[k]==kk)])}
      SNP_table[i,"AAfrom"] = ko
      ko=c()
      kk = SNP_table[i,"AAto_codon"]
      for (k in 1:6) {ko=c(ko,rownames(AAtrans)[which(AAtrans[k]==kk)])}
      SNP_table[i,"AAto"] = ko
    }
  }
  SNP_table[!is.na(SNP_table$protein) & SNP_table$AAfrom == SNP_table$AAto,"syn"] = "Syn"
  SNP_table[!is.na(SNP_table$protein) & SNP_table$AAfrom != SNP_table$AAto,"syn"] = "NoSyn"
  SNP_table[which(SNP_table$syn == "NoSyn"),"mutation"] = paste0(SNP_table[which(SNP_table$syn == "NoSyn"),"AAfrom"],SNP_table[which(SNP_table$syn == "NoSyn"),"AAloc"],SNP_table[which(SNP_table$syn == "NoSyn"),"AAto"])
  SNP_table[which(SNP_table$syn == "NoSyn"),"mutation"] = paste0(SNP_table[which(SNP_table$syn == "NoSyn"),"protein"],":",SNP_table[which(SNP_table$syn == "NoSyn"),"mutation"])
  SNP_table = SNP_table[c(`SNP`,c("protein","AAfrom","AAloc","AAto","syn","mutation"))]
  return(SNP_table)
}

