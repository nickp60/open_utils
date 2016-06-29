# Illumina summary csv to Metphylian
#  Nicholas Waters 20160418
################  gather inputs, read in, and aggregate into a single list
DEBUG=T
if(DEBUG){
  args<-c("~/SJ_data/Lakes/","Aggregate Summary-31540571/") 
#  args<-c("~/SJ_data/Anacostia/")
} else{
  args<-commandArgs(T)
}
if(!length(args)>=1){
  print("usage: Rscript illumina2metphylian.R path/to/input/directory/")
  stop("see usage, try again!")
}
input_directory<-args[1]
dirs_to_exclude<-ifelse(is.na(args[1]), "standin", c(args[2:length(args)]))
all_files<-
  grep(".*summary.csv", #get summary files
       grep(paste(dirs_to_exclude,".*",sep=""), #except ones in excluded paths
            list.files(input_directory, recursive = T), #from the input directory
            invert=T, value=T)
       , value=T)
aggregated<-lapply(all_files, function(x){
  c<-read.csv(paste(input_directory,x,sep=''), stringsAsFactors = F)
  c$sample<- gsub("(.*)\\/(.*)\\.summary\\.csv", "\\2", x)
  return(c)
})
agg_df<-do.call(rbind, aggregated)
################   clean up column names
colnames(agg_df)[grep("X\\._*", colnames(agg_df))]<-"percent_hits"
################  rename, reduce, remove raw hit number
agg_df$ID<-paste("k__",agg_df[,1], 
                 "|p__",agg_df[,2], 
                 "|c__",agg_df[,3], 
                 "|o__",agg_df[,4], 
                 "|f__",agg_df[,5], 
                 "|g__",agg_df[,6], 
                 "|s__",agg_df[,7], sep="")
agg_df<-agg_df[,8:ncol(agg_df)]
agg_df[grep("num_hit", colnames(agg_df))]<-NULL # remove raw hit, retain percentages
################   begin conversion to format in "mergend_abundance_table.txt
# with reshape
if(any(grepl("reshape2", installed.packages()))){
  library(reshape2)
  cast_aggregated<-dcast(agg_df, ID~sample, value.var = "percent_hits", fill = 0)
} else{
  cast_aggregated<-dcast(agg_df, ID~sample, value.var = "percent_hits", fill = 0)
  tempdf<-agg_df[agg_df$sample==unique(agg_df$sample)[1],]
  tempdf$sample<-NULL
  tempdf$percent_hits<-NULL
  for (i in unique(agg_df$sample)){
    tempdf<-merge(tempdf, agg_df[agg_df$sample==i,c("ID","percent_hits")], by="ID", all=T)
    colnames(tempdf)[grep("percent_hits", colnames(tempdf))]<-i
  }
  tempdf[is.na(tempdf)]<-0
  cast_aggregated<-tempdf
}
### writ eout the results
write.table(cast_aggregated, 
            file = paste(input_directory,"formated_aggregated_table.txt",sep=''),
            row.names=FALSE,sep = "\t")




