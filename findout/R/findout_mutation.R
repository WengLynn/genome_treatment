
findout_mutation <- function(filename,key_pattern ) {
  # coln -------------------------------------------------------------------
  library(tidyr)
  library(stringr)
  ref_file = "reffa.fasta"
  cod = paste0("grep -A1 '",key_pattern,"' ",filename," > ",ref_file)
  system(cod)
  reffa =read.table(ref_file,sep = "\n")
  reffa = as.character(reffa[2,1])
  mmref = unlist(strsplit(reffa,""))
  mmref_withoutgap = mmref[!mmref=="-"]
  withoutgap = which(mmref=="-")
  mmcount = length(mmref)
  ref=ifelse(mmref=="-",1,0)
  barcount = sum(ref)
  insertpoi = which(ref==1)
  coln = c(1:(mmcount-barcount))
  for (i in 1:barcount) {
    p = paste((insertpoi[i]-i),"_",(insertpoi[i]-i+1),"insNo.",i,sep = "")
    coln = append(coln,p,(insertpoi[i]-1))
  }
  loc = coln
  loc[1:length(mmref)] = 1
  loc[which(mmref == "-")]=0
  loc = cumsum(loc)
  loc[which(mmref == "-")]=0
  # loop --------------------------------------------------------------------
  write.table(t(as.data.frame(c("Type","loc_from","loc_to","NNfrom","NNto","seqname"))),"indel_info.txt",col.names = F,row.names = F,sep = "\t",quote = F,append = T)
  write.table(t(as.data.frame(c("Type","loc_from","loc_to","NNfrom","NNto","seqname"))),"SNP_info.txt",col.names = F,row.names = F,sep = "\t",quote = F,append = T)
  lineCnt = 0
  con <- file(filename, "r")
  while (1) {
    oneline = readLines(con, n = 1)
    if(length(oneline) ==0){break}
    if(length(grep(">",oneline)) == 1){name =oneline;namestr = unlist(str_split(name,"\\|"))}#; print(name)} #; name =unlist(strsplit(name,"\\|"))}
    if(length(grep(">",oneline)) == 0) {
      #   write.table(paste(name,lineCnt,sep = "\t"),"namecheck.txt",sep = "\n",quote = F,col.names = F,row.names = F,na = "",append = T)
      mm = unlist(strsplit(oneline,""))
      mmcount = ifelse(mm == "N" ,1,0)
      mmcount2 = ifelse(mm %in% c("-","A","T","C","G","N") ,0,1)
      mmsum = sum(mmcount)
      mmsum2 = sum(mmcount2)
      if(mmsum<15 & mmsum2<50){
        write.table(lineCnt,"QC_line.txt",sep = "\n",quote = F,col.names = F,row.names = F,na = "",append = T)
        lineCnt = lineCnt+1
        if (lineCnt%%10000==0){print(lineCnt)}
        mm2=unlist(strsplit(as.character(oneline),""))
        mmc =mm2
        mmref_temp=mmref
        if(length(mm)!=length(mmref)){mmref_temp = mmref[1:length(mm)]}
        mmref_temp[!mm2 %in% c("-","A","T","C","G")] = mm2[!mm2 %in% c("-","A","T","C","G")]
        mmc[mm2 == mmref_temp]= 100                               #same = 100
        mmc[mm2 != mmref_temp & mm2 =="-"]= 5                     #dele = 5
        mmc[mm2 != mmref_temp & mmref_temp =="-"]= 1                   #inse = 1
        mmc[mm2 != mmref_temp & mm2 != "-" & mmref_temp != "-"] = 10  #SNP = 10
        mm_del = which(mmc == 5)
        mm_ins = which(mmc == 1)
        mm_SNP = which(mmc == 10)

# deletion ----------------------------------------------------------------


        if(length(mm_del)>0){
          indel_table=""
          mm_out = mm2[-withoutgap]
          loc_del = loc[mm_del]
          loc_len = length(loc_del)
          if (length(mm_del)!=1){
          loc_s = c(loc_del[1],loc_del[1+which(loc_del[2:loc_len]-loc_del[1:(loc_len-1)]!=1)])
          loc_e = c(loc_del[which(loc_del[2:loc_len]-loc_del[1:(loc_len-1)]!=1)],loc_del[loc_len])
          NNf  = NNe = c()
          for (del_ty in 1:length(loc_s)) {
            NNf =c(NNf, paste(mmref_withoutgap[loc_s[del_ty]:loc_e[del_ty]],collapse = ""))
            NNe =c(NNe, paste(mm_out[loc_s[del_ty]:loc_e[del_ty]],collapse = ""))
          }
          indel_table = cbind("del",loc_s,loc_e,NNf,NNe,name)
          write.table(indel_table,"indel_info.txt",col.names = F,row.names = F,sep = "\t",quote = F,append = T)
          loc_del=0
          }
          if(length(mm_del)==1){
            loc_s = loc_del
            loc_e = loc_del
            NNf  = NNe = c()
            NNf =mmref_withoutgap[loc_del]
            NNe =mm_out[loc_del]
            indel_table = cbind("del",loc_s,loc_e,NNf,NNe,name)
            write.table(indel_table,"indel_info.txt",col.names = F,row.names = F,sep = "\t",quote = F,append = T)
            loc_del=0
          }

        }

# insertion  --------------------------------------------------------------

        if(length(mm_ins)>0){
          indel_table=""
          mm_out = mm2[-withoutgap]
          loc_del = loc[mm_ins]
          loc_len = length(loc_del)
          if (length(mm_ins)!=1){
            loc_s = c(loc_del[1],loc_del[1+which(loc_del[2:loc_len]-loc_del[1:(loc_len-1)]!=1)])
            loc_e = c(loc_del[which(loc_del[2:loc_len]-loc_del[1:(loc_len-1)]!=1)],loc_del[loc_len])
            NNf  = NNe = c()
            for (del_ty in 1:length(loc_s)) {
              NNf =c(NNf, paste(mmref_withoutgap[loc_s[del_ty]:loc_e[del_ty]],collapse = ""))
              NNe =c(NNe, paste(mm_out[loc_s[del_ty]:loc_e[del_ty]],collapse = ""))
            }
            inindel_table = cbind("ins",loc_s,loc_e,NNf,NNe,name)
            write.table(indel_table,"indel_info.txt",col.names = F,row.names = F,sep = "\t",quote = F,append = T)
            loc_del=0
          }
          if(length(mm_ins)==1){
            loc_s = loc_del
            loc_e = loc_del
            NNf  = NNe = c()
            NNf =mmref_withoutgap[loc_del]
            NNe =mm_out[loc_del]
            indel_table = cbind("ins",loc_s,loc_e,NNf,NNe,name)
            write.table(indel_table,"indel_info.txt",col.names = F,row.names = F,sep = "\t",quote = F,append = T)
            loc_del=0
          }

        }

# SNP ---------------------------------------------------------------------


        if(length(mm_SNP)>0){
          write.table(cbind("SNP",loc[mm_SNP],mmref[mm_SNP],mm2[mm_SNP],name),"SNP_info.txt",col.names = F,row.names = F,sep = "\t",quote = F,append = T)
        }
      }
    }
  }
  close(con)



}
