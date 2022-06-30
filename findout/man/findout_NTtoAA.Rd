\name{findout_NTtoAA}
\alias{findout_NTtoAA}
\title{findout_NTtoAA}
\usage{
findout_NTtoAA(Seq_fasta ,ORF_table_file,SNP_table_file,ORFStartCol = "Start",ORFEndCol = "End",ORFGeneCol="Gene",SNP = "SNP")
}
\description{
Transfer nucleotide mutations to Amino acid mutations.
}
\examples{
findout_NTtoAA("reference.fasta","ORF.csv","SNP.csv")
}
