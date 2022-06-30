shorten_site <- function(filename,key_pattern,Start=300,End=29000,hm =0) {
  # coln -------------------------------------------------------------------
  library(tidyr)
  library(stringr)
  ref_file = "reffa.fasta"
  cod = paste0("grep -A1 '",key_pattern,"' ",filename," > ",ref_file)
  system(cod)
  reffa =read.table(ref_file,sep = "\n")
  reffa = as.character(reffa[2,1])
  mmref = unlist(strsplit(reffa,""))
  mmcount = length(mmref)
  ref=ifelse(mmref=="-",1,0)
  barcount = sum(ref)
  insertpoi = which(ref==1)
  coln = c(1:(mmcount-barcount))
  for (i in 1:barcount) {
    p = paste((insertpoi[i]-i),"_",(insertpoi[i]-i+1),"insNo.",i,sep = "")
    coln = append(coln,p,(insertpoi[i]-1))
  }

  # loop 1--------------------------------------------------------------------
  con <- file(filename, "r")
  cnt_seq  = 0
  cnt_seq[1:mmcount] = 0
  while (1) {
    oneline = readLines(con, n = 1)
    if(length(oneline) ==0){break}
    if(length(grep(">",oneline)) == 0) {
      mm = unlist(strsplit(oneline,""))
      mm[which(mm == "N"|mm=="n")]<- mmref[which(mm == "N"|mm=="n")]
      cnt_seq[which(mm!=mmref)]=cnt_seq[which(mm!=mmref)]+1
    }
  }
  close(con)
  colnstart = which(coln == Start)
  colnend =which(coln == End)
  coln = coln[colnstart:colnend]
  cnt_seq = cnt_seq[colnstart:colnend]
  coln =coln[which(cnt_seq>hm)]
  write.table(coln,paste0(filename,"_mic_loc.txt"),col.names = F,row.names = F,sep = "",quote = F,append = F)

# loop2 -------------------------------------------------------------------



  con <- file(filename, "r")
  while (1) {
    oneline = readLines(con, n = 1)
    if(length(oneline) ==0){break}
    if(length(grep(">",oneline)) == 1){write.table(oneline,paste0(filename,"_",Start,"_",End,"_mic.fasta"),col.names = F,row.names = F,sep = "\t",quote = F,append = T)}
    if(length(grep(">",oneline)) == 0) {
      mm = unlist(strsplit(oneline,""))
      mm = mm[colnstart:colnend]
      mm= mm[which(cnt_seq>hm)]
      mm =paste(mm,collapse = "")
      write.table(mm,paste0(filename,"_",Start,"_",End,"_mic.fasta"),col.names = F,row.names = F,sep = "",quote = F,append = T)
    }
  }
  close(con)
}


