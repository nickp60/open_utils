#! /usr/bin/Rscript
#############################################################################################
DEBUG=F
source("~/GitHub/open_utils/NargParse/NargParse.R")
#############################################################################################
help<-"
Blast Virulence Parser
20160616

Minor version changes:
 - PCA
 - Plasmid mapping

USAGE:
REQUIRED ARGUMENTS:
[-r resultsfromdcmegablast.csv ] 
[-o output_firectory]	    path to output directory. 
OPTIONAL ARGUMENTS:
[-t <1-100>] bit score percentage similarity to call a hit
[-p genome_plasmid_accessions.txt] file used with get_genomes.py
[-verb T or F]	          Turn verbose output on and off


"
require(dplyr)
require(reshape2)
require(ggplot2)
this_version<-'0.5'
test_args<-c(
  "-r", "~/GitHub/FB/Ecoli_comparative_genomics/results/20161122_snagnblast/snagnblast_results_with_ref/20161122_results_merged.csv", # path to blast results
  "-t","25", # threshold percentage
  "-p","~/GitHub/FB/Ecoli_comparative_genomics/data/pathogenic_e_coli_genome_accessions.txt", #for genomes with plasmids
  "-o","~/GitHub/FB/Ecoli_comparative_genomics/results/blast_virulence_parser_output/20161122/",  # output directory for figs etc
  "-verb", "T") #Verbose Bool
if (DEBUG | length(commandArgs(T))==0){
  print("Caution! Using with debug dataset")
  pre_args<-test_args
} else{
  pre_args<-commandArgs(T)
}
## @knitr part1

# parse the commandline arguments 
cline_args<-parse_flagged_args(x = pre_args, required = c("-r", "-o"), 
                               optional = c("-t", "-verb", "-p"),
                               version = this_version, help_message = help, test_args=test_args,DEBUG=F)
VERBOSE_bool<-ifelse(arg_bool(cline_args['verb']), T,F)

blast_results<- arg_file(cline_args['r'])

#genomes_with_plasmids<-arg_file(cline_args["p"])

perc_similarity_thresh<-arg_numeric(cline_args["t"], default = 50)*.01
if(perc_similarity_thresh>1 |perc_similarity_thresh<=0){stop("percentage must be between 1 and 100")}
datetime<-gsub("[:-]","",gsub("([\\d:-]*)(\\D*)","\\1",Sys.time()))

output_directory<-arg_directory(cline_args["o"], makedir = T, need_write = T, recursive = F)
output_pre<-paste0(datetime,"_blast_virulence_parser_output_")
setwd(output_directory)
##########
#read in csv
data<-read.csv2(blast_results, stringsAsFactors = F, sep=",")

#make bit_score numeric
data$bit_score<-as.numeric(data$bit_score)
# make column for ifelse  value/maxvalue>threshold, or if it is threshold or greater than the max bit score per gene.
data2<-data %>%
  group_by(query_id) %>%
  mutate(norm_bit_score=(bit_score/max(bit_score))*100)%>%
  mutate(pass=ifelse((bit_score/max(bit_score))>perc_similarity_thresh, 1,0))%>%
  filter(pass==1)

#simplify for casting to long-format data
cast_data = data.frame(gene=paste0(data2$accession,"_", data2$name), genome=data2$subject_id, pass=data2$pass, norm_bit_score=data2$norm_bit_score, stringsAsFactors = F)

# combine contigs
# TODO if this could be fixed prior it would be way less fragile
cast_data$genomesimple<-
  gsub("(.*\\d)(_c.*)", "\\1", gsub("(.*)ref\\|(.*)\\|", "\\2", cast_data$genome))
cast_data$genome<-NULL


if(!is.na(cline_args["p"])){  #if given a plasmid file
  gwp<-readLines(cline_args["p"]) # genomes_with_plasmids) # read it in
  #remove doubled hashes, so that only the one preceding the aaccesion will be used for naming
  for(i in 1:length(gwp)){
    if(grepl("^#.*", gwp[i]) && grepl("^#.*", gwp[i+1])){
      gwp[i]<-NA
    }
  }
  renaming_df<-data.frame(genomesimple=NA, genome=NA, stringsAsFactors = F) #creat recipient DF
  renaming_index<-1 #counter
  gwp<-gwp[complete.cases(gwp)] #get rid of extra hash rows marked above
  gwp[length(gwp)+1]<-"#" #add buffer hash to end of 
  hashes<-grep("^#.*", gwp) #indexes of all hashes
  for(i in 1:(length(hashes)-1)){#print(i)   # for each hash (excluding the dummy) index,
    for(j in (hashes[i]+1):(hashes[i+1]-1)){  # and for each line inbetween the hash and the next,
      new_line<-data.frame(genome=gsub(" ","_", gsub("#","", gwp[hashes[i]])), #make a datataframe with the genome name and the included accession
                           genomesimple=gwp[j], stringsAsFactors = F)
      renaming_df<-rbind(renaming_df, new_line) #add it to the recipient structure
    }
  }
  renaming_df<-renaming_df[complete.cases(renaming_df),] #remove first NA row
  cast_data<-merge(cast_data, renaming_df, by="genomesimple", all.x=T) #add genome column to cast_data
  cast_data$genome<-ifelse(is.na(cast_data$genome), cast_data$genomesimple, cast_data$genome) #tidy up names by copying over the non-reannotated names from genomesimplify column
  cast_data$genomesimple<-NULL #no longer needed
  control_strains<-renaming_df$genome #this is used when selecting only genes returned in at least one of the positive control strains
  
} else{
  control_strains<-unique(grep("NC_*", cast_data$genomesimple, value=T))  #control strains if just using the NC_* genome names
  renaming_df<-control_strains
  cast_data$genome <- cast_data$genomesimple
  
}

  
##########
#remove duplicates (usually hits on multiple redundant contigs)
cast_data<-cast_data[!duplicated(cast_data),c("gene", "genome", "pass", "norm_bit_score")]
#get best hit per genome, selecting only those with hits
cast_data2<-cast_data%>%
  group_by(gene, genome)%>%
  mutate(best=ifelse(norm_bit_score==max(norm_bit_score),"best","notbest"))%>%
  filter(best=="best")
cast_data2$best<-NULL #remove max column
# make wide again for exporting table, etc
wide_results<-dcast(data=cast_data2[,colnames(cast_data2) != "pass"], 
                   formula = gene~genome, 
                   value.var = "norm_bit_score", fill = 0) #, fun.aggregate = length)
row.names(wide_results)<-wide_results$gene
wide_results$gene<-NULL

#wide_results<-  as.data.frame(ifelse(wide_results>=1, 1,0), stringsAsFactors = F) #ignore duplicates; yes, I hate myself
wide_results$gene<-row.names(wide_results)

#  remove genes not found in any control strains; 
#head(wide_results)
#colnames(wide_results)
controls<-c("gene", control_strains)
genes_ignored<-c()
for (i in wide_results$gene){
  if(sum(wide_results[wide_results$gene==i,control_strains])==0 ){
    print(paste("removing all",i, "from dataset" ))
    genes_ignored<-c(genes_ignored,i)
    wide_results<-wide_results[wide_results$gene !=i , ]
  }
}
#TODO write out to csv for supplementary data?
write.csv(wide_results, file = paste0(output_pre,"hit_table.csv"))
cat(
  (c("##########  OUTPUT from blast_virulence_parser.R   ########",
    this_version,
    "Command Line Arguments used:",
    cline_args,
    "Control genomes and their plasmids:",
    (renaming_df$genomesimple),"",(renaming_df$genomes),
    "Genes ignored for analysis:",
    genes_ignored)),file = paste0(output_pre,"runinfo.txt"), sep = "\n", append = F )

#add clustering
set.seed(27)

require(gplots)
colors = c(seq(-1,49,length=10),seq(50,100,length=10))
my_palette <- colorRampPalette(c("#FFEBEB", "#FF2B2B"))(n = 19)
pdf(file = paste0(output_pre,"heatmap.pdf"), width = 10, height = 14)
par(cex.main=0.6)
heatmap.2(t(as.matrix(wide_results[,!colnames(wide_results) %in% c("k_clus", "gene")])),
          col=my_palette, key.ylab = "Count", key.xlab="Normalized Bit Score",
          breaks=colors, density.info='histogram', tracecol = 'black',
          Rowv=T, Colv=T, 
          dendrogram="both", densadj=.5, symkey=F, #col=redblue(15), 
          key=T, keysize=1,
          trace="none",  
          offsetRow=-.3, offsetCol=-.3,
          cexRow  = 0.5, cexCol = 0.5,
          adjRow=c(0,0), adjCol=c(1,.5),
          margins=c(6.5,6.5)
)
dev.off()









