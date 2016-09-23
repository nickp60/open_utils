#  Semi-globl Alignment vis
args<-commandArgs(T)
if(length(args) != 2){
  print("Usage: Rscript gapped_alignment_vis.R file1.fasta file2.fasta")
  stop("Must have two and only two arguments")
}
suppressPackageStartupMessages(require("Biostrings"))
pattern <- readDNAStringSet(args[1], format="fasta")
subject <- readDNAStringSet(args[2], format="fasta")
#subject <- readDNAStringSet("~/riboSeed_sandbox/seed2/results/SPAdes_20160815_ribosomal_RNA_16S226621_228162_riboSnag/contigs.fasta", format="fasta")
#pattern <- readDNAStringSet("~/GitHub/FB/Ecoli_comparative_genomics/scripts/riboSeed_pipeline/batch_coli/20160815_ribosomal_RNA_16S226621_228162_riboSnag.fasta", format="fasta")
# subject <- readDNAStringSet("~/Desktop/tempgenome.fasta", format="fasta")
#pattern <- DNAString("ATGCCCATATA")
#subject <- DNAString("ATATGGGCAT")
# 

print_alignment<-function(subject, pattern, breakp=60){
  # switch pattern and subject if pattern is longer
  if(nchar(pattern) > nchar(subject)){
    print("switching subject and query as query was longer")
    temp <- pattern
    pattern <- subject
    subject <- temp
  }
  #  sort so that start(subject(x)) >1
  x<-(pairwiseAlignment(pattern = pattern, subject = subject, type="overlap"))
  y <- (pairwiseAlignment(pattern = reverseComplement(pattern), subject = subject, type="overlap"))
  if (score(x)>=score(y)){
    x <- x
  } else{
    print("using the reverse compliment")
    x <- y
    pattern<-reverseComplement(pattern)
  }
  if (score(x) < (.1*nchar(pattern))){
    print("No match found with a score above 1/10th the length of the pattern")
    return(1)
  } else{
    print(score(x))
  }
  if (start(subject(x))==1){
    print("switching pattern and subject for convenience")
    temp=subject
    subject=pattern
    pattern=temp
    x<-(pairwiseAlignment(pattern = pattern, subject = subject, type="overlap"))
  }
  max<-ifelse(nchar(pattern(x))>nchar(subject(x)), nchar(pattern), nchar(subject))
  # print(as.character(x))
  new_query<-paste0(paste0(rep("-", start(subject(x))), collapse = ""), 
                    pattern(x),
                    substring(pattern, end(pattern(x)), nchar(pattern)),
                    paste0(rep("-", max-nchar(as.character(x))), collapse = ""))
  new_subject<-paste0(paste0(rep("-", start(pattern(x))-1), collapse = ""),
                      substring(subject,1,start(subject(x))),
                      subject(x), 
                      paste0(rep("-", (nchar(pattern)-end(pattern(x)))), collapse = ""))
  consensus = ""
  for (j in 1:nchar(new_query)){
    # que<-substring(new_subject, (i*breakp)+1, (breakp*i)+breakp)[j]
    # subsubstring(new_query,   (i*breakp)+1, (breakp*i)+breakp)[j]
    if (substring(new_subject,j ,j) == substring(new_query,j ,j)){
      consensus<-paste0(consensus, "|")
    } else{      
      consensus<-paste0(consensus, " ")
    }
  }
  for (i in 0:(round(max/breakp))){
    cat(paste0(substring(new_subject, (i*breakp)+1, (breakp*i)+breakp), "\n"))
    cat(paste0(substring(consensus, (i*breakp)+1, (breakp*i)+breakp), "\n"))
    cat(paste0(substring(new_query, (i*breakp)+1, (breakp*i)+breakp), "\n\n"))
  }
  return(0)
}

for (i in 1:length(subject)){
  print(names(subject[i]))
  print_alignment(subject=subject[i], pattern=pattern, breakp = 80)
}
#

