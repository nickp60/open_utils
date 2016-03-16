#!/usr/bin/Rscript
DEBUG<-F
################################################################################
################################################################################
#                           gbparseR
#
#                         by Nick Waters
#                           20160325
#                         Version 0.3.93
################################################################################
################################################################################

#   Usage: $ Rscript gbparseR.R input.gb  output_path *upstream_ and_downstream_bps
#                                                     *optional

# Minor update 0.3.93 
#fixed NA's in strand

#  known bugs:  need to make support for joined CDS's , ie:
# "complement(join(332323..2323,223234.2342325.))
#  -broken gff3 output.  THIS IS THE BIGGEST TODO


# Cant have upstream region definitions and no output path.  you will break it. :(

#  the functions will extract 
# informaton from .gb genbank files, uing a lot of regular expressions and the 
# readlines function.  Everything is base R, so that should make your life easier 
# when trying to manage packages..
#
# get sequences from NCBI ising the dropdown menu for exporting a sequence 
# after ensuring that you are displaying the full sequence at the bottom of the page.

#DANGER:  previous versions have a 7bp frame deletion (cause I'm a silly and make silly mistakes)
################################################################################
################################################################################
#Test files
#source.file = "~/GitHub/BlastDBs/CP000253_8325.gb"
#source.file = "~/GitHub/BlastDBs/N315.gb"
#source.file = "~/GitHub/BlastDBs/uams1.gb"
#source.file ="~/GitHub/BlastDBs/FPR3757_LAC.gb"
#source.file ="~/GitHub/BlastDBs/TCH1516.gb"
if(DEBUG){
  args<-c("~/uams1_rs.gb","../gbparse_output/", "500")
} else{
  args<-commandArgs(TRUE)
}
if (is.na(args[1]) | !file.exists(args[1])){  #if input is blank or the file doesnt exist,
  stop("Error: needs an input file!")
}
source.file<-args[1]
if (is.na(args[2])){
  print("using current working directory and default region size of 500")
  cur_dir<-getwd()
  working_dir<-paste(cur_dir,"/", gsub("-","",Sys.Date()),"gbparseR","/", sep="")
  if (!dir.exists(working_dir)){dir.create(working_dir)}
  setwd(working_dir)
} else if (is.numeric(args[2])){
    stop("cannot give region width without an output directory")
} else {
    output_path<-args[2]
    if (!dir.exists(output_path)){dir.create(output_path)}
      setwd(output_path)
      working_dir<-output_path
}
if (is.na(args[3]) | is.na(as.numeric(args[3])) | 0>(as.numeric(args[3]))){
  print( "no region or non-numeric region; using 500bp as region width")
  upstream<-downstream<-500
} else{
  downstream<-upstream<-as.numeric(args[3])
}
#downstream<-upstream #same upstream as downstream.  need more flexibility? figure it out..
#  sanity check:
print(paste("genome source file: ", source.file))
print(paste("output directory: ", working_dir))
print(paste("upsteam and downstream region width: ", upstream))
################################################################################
####  Extract each scaffold and write out to files in a new 
####  directory
split_if_scaffolded<-function(source.file){
  input <- readLines(source.file)
  if(length(grep("LOCUS", input))>1){
    print("splitting scaffolds into two files for easier handling...") 
#    scaffs<-length(grep("LOCUS", input))
    scaffCoords<-c(grep("LOCUS", input), length(input))
    for( i in 1:(length(scaffCoords)-1)){#print(i)}
      scaf<-
        gsub(" *","",
             gsub("VERSION (.*) (.*)","\\1", grep("VERSION", input, value = T)[i]))
      scaffFileName<-paste(scaf,"_scaf_", i, ".txt", sep="")
      ifelse (i==1,
              splitInput<-input[scaffCoords[i]:scaffCoords[i+1]-1],
              splitInput<-input[scaffCoords[i]:scaffCoords[i+1]])
      print(scaffFileName)
      write(splitInput, scaffFileName)
    }
  } else{
    scaf<-gsub(" *","", 
               gsub("VERSION (.*) (.*)","\\1", grep("VERSION", input, value = T)))
    print(paste("saving formatted scaffold:",paste(scaf, "_scaf_1.txt", sep="")))
    write(input, paste(scaf, "_scaf_1.txt", sep=""))
  }
  gc()
}
#^^^^^^^^^  
#split_if_scaffolded(source.file)

################################################################################
####  extract metadata and meat ( the actual features).  spits out a list 
#     containing:
#         the meta info
#         The features
#         the nucleotide sequence
load_gb<-function(source.file){
  input <- readLines(source.file)
    meta<-input[grep("LOCUS", input)[1]:(grep("FEATURES", input)[1]-1)]
    meat<-input[!input %in% meta]
    meat<-meat[grep("FEATURES", meat):length(meat)]
    ifelse(any(grepl("ORIGIN", meat)) ==F | length(grep("ORIGIN", meat))>1, 
           "Error!  You have no ORIGIN, or more than one ORIGIN\n
                 Do you have the complete sequence?","" )
    seq<-meat[grep("ORIGIN", meat): length(meat)]
    meat<-meat[grep("FEATURES", meat): grep("ORIGIN", meat)] #get rid of sequence from features
    return(list(meta, meat, seq))
#  }
}
#^^^^^^^^^  
#returns<-load_gb(source.file)
################################################################################
#####  clean up sequence
clean_sequence<-function(seq){#} (seq in grep("seq\\d", ls(), value=T)){#print(seq)}
  fullraw<-seq #get(seq)
  nospace<-gsub("\\s", "",fullraw[2:length(fullraw)])
  nonum<-paste(gsub("\\d", "",nospace), sep="", collapse="")
  nonum<-gsub("(.*)(//.*)","\\1",nonum)
  nonum
}

#^^^^^^^^^ 
#returns[[3]]<-clean_sequence(seq = returns[[3]])

################################################################################
###  start making  dataframe containing the location coordinates
get_ranges<-function(x){
  if (grep("FEATURES", x)!=1){stop}   #  make sure its 
  x<-c(x, " 27..27 ")  #  hey there, buffer row
  # extract ranges 
  featList<-  #  note!  this removes ">" exception from locus; see exception in tag"
    gsub("\\)|\\,|\\(|-","", 
         gsub("(\\D*)(\\d+\\.{2}[\\>,0-9]{1}\\d*)(\\)*.*)", "\\2", 
              
              grep("(\\d*\\.{2}\\d*)",x, value=T)))
  #featList<-gsub("#clean ranges
  indexList<-grep("(\\d+\\.{2}[\\>,0-9]{1}\\d*)",x)
  preZ<-data.frame(index= indexList, loc=featList, stringsAsFactors = F)
  # this next bit gets rid of duplicates;  neat, huh?
  for (i in 2:length(preZ$loc)){
    if (preZ[i-1,"loc"]==preZ[i, "loc"]){preZ[i,"loc"]<-"NA"}
  }
  preZ<-preZ[preZ$loc != "NA",]
  z<-preZ   
  z$id<-1:nrow(z) #relevel id's
  z$loc_start<-as.numeric(gsub("(\\d*)(\\.\\.[\\>,0-9]{1}\\d*)", "\\1", z$loc))
  z$loc_end<-as.numeric(gsub (">","", gsub("(\\d*\\.\\.)(\\d*)", "\\2", z$loc)))
  z$next_loc<-unlist(c(lapply(z$id[1:nrow(z)-1], function(w){       #  get start of next locus
    as.numeric(z[z$id==w+1, "loc_start"])}), 0))
  z$loc<-NULL
  z
  
}
#^^^^^^^^^ 
#ranges<-get_ranges(returns[[2]])


################################################################################
#  this removes references to locations from the locations lis, and 
#  relevels the id's to get continuous numbering

clean_ranges_rna<-function(gl, data){
  for (i in gl[1,"id"]:gl[nrow(gl) - 1, "id"]){
    j<-i+1
    entry<-
      data[gl[gl$id==i,"index"]:gl[gl$id==j,"index"]]
    gl[gl$id==i, "keep"]<-ifelse(
      grepl("anticodon=|transl_except=", entry[1]), "NA", "keep"
    )
  }
  gl[gl$id==max(gl$id, na.rm = T), "keep"]<-"keep"    #  take care of the buffer row;  need this to be included
  gl<-gl[gl$keep != "NA",]
  gl[,"keep"]<-NULL
  gl[,"id"]<-1:nrow(gl)#relevel id to get consec numbers
  gl
}

#^^^^^^^^^ 
#ranges<-clean_ranges_rna(ranges, returns[[2]])

################################################################################

#  this extracts the features, augmenting the ranges dataframe with gene info, 
#  type, etc
extract_features<-function(data, ranges, locus_tag="locus_tag", debug=F){ #y=grep("\\sgene\\s", x)[1]
  
  #   check for buffer row, stop if  absent
  if(tail(ranges,1)$loc_start!=27 | tail(ranges,1)$loc_end != 27){
    stop("where's the buffer row? did you forget to get and clean the ranges?")
  }
  z<-ranges
  for (i in ranges[1,"id"]:(ranges[nrow(ranges),"id"]-1)){ #  print(i)}
    if(debug){print(i)} # print lines
    #  can this be made into an apply?
    j<-i+1
    hits<-grep("\\d+\\.\\.[\\>,0-9]{1}\\d*",  # ">" is part of a new convention for stuff?
               data[ranges[ranges$id==i,"index"]:ranges[ranges$id==j,"index"]])
    #  this defines the index of where to look for the entry header;  default is 1:3 
    #  for non-gene features;  for genes, it finds the loc pattern \\d..\\d in the 
    #  region and bound of the entry as the second, middle occurance.
    entry_end <-ifelse(length(hits)>=3,
                       hits[2], 4)   ####################   observe how many lines 
                                     #################### are b/w "gene" and "CDS"
    entry<-data[ranges[ranges$id==i,"index"]:ranges[ranges$id==j,"index"]]
    entry_head<-entry[1:entry_end]
    #}#  quick check
    if(!any(grepl("source|repeat_region|STS|tmRNA|assembly_gap|misc_feature|misc_binding|CDS|misc_RNA|rRNA|ncRNA|tRNA", 
                  entry_head))){
      if(any(grepl("pseudo", entry))){
        print(paste("caution! id =", i, "is possibly a pseudogene"))
      } else {
        print(paste(" uh oh...  we got a rogue entry:",i,". Skipping..."))
        next()
      }
      #
    }
    # go time;  stuff is arranged most frequent to least.  probs will need
    #to add lotsa "try" stuff to this
    if(any(grepl("CDS", entry_head))){
      z[z$id==i, "type"]<-"CDS"
      
      pre_locus<-gsub(paste("(.*\\s*\\",locus_tag,"=)(.*)(\")", sep=""),"\\2", 
                      grep(locus_tag,entry, value=T )[1])
      try(z[z$id==i, "gene"]<-
            gsub("(.*\\s*\\gene=\")(.*)(\")","\\2", 
                 grep("gene=",entry, value=T )[1]), silent=T)
      try(z[z$id==i, "old_locus_tag"]<-
            gsub("(.*\\s*\\old_locus_tag=\")(.*)(\")","\\2", 
                 grep("old_locus_tag=",entry, value=T )[1]), silent=T)
      
      z[z$id==i,"locus_tag"]<-
        gsub("(.*?\")(.*)", "\\2",  pre_locus)
      try(z[z$id==i, "inference"]<-
            gsub("(.*\\s*\\inference=)(.*\")(.*)","\\3", 
                 grep("inference",entry, value=T )),silent=T)
      #    z[i, "sequence_ref"]<-
      #      gsub("(.*\\s*\\sequence:)(.*)","\\2", grep("sequence:",entry, value=T ))
      z[z$id==i, "direction"]<-ifelse(!any(grepl("complement",entry_head)), "leading", "compliment")
      try(z[z$id==i, "note"]<-
            gsub("  ","",paste0(entry[grep("note=", entry): 
                                        (grep("codon_start",entry)-1)], collapse="  ")),
          silent=T)
      try(z[z$id==i, "codon_start"]<-
            gsub("(.*\\s*\\codon_start=)(.*)","\\2", 
                 grep("codon_start",entry, value=T )), silent=T)
      try(z[z$id==i, "transl_table"]<-
            gsub("(.*\\s*/transl_table=)(.*)","\\2", 
                 grep("transl_table",entry, value=T )), silent=T)
      try(z[z$id==i, "product"]<-
            gsub("(.*\\s*/product=)(\")(.*)(\")","\\3", grep("product",entry, value=T )),
          silent=T)
      try(z[z$id==i, "protein_id"]<-
            gsub("(.*\\s*/protein_id=)(\")(.*)(\")","\\3", 
                 grep("protein_id",entry, value=T )), silent=T)
      try(z[z$id==i, "db_xref"]<-
            gsub("(.*\\s*/db_xref=)(\")(.*)(\")","\\3", 
                 grep("db_xref",entry, value=T )), silent=T)
      try(z[z$id==i,"translation"]<-
            paste(gsub("([^A-Z]*)([A-Z]*)([^A-Z]*)","\\2",
                       gsub("translation=", "", entry[grep("translation=", entry):length(entry)])), 
                  collapse="",sep=""), silent=T)
      #####  
    } else if(any(grepl("pseudo", entry[length(entry)-1]))){
      z[z$id==i, "type"]<-"pseudo"
      pre_locus<-gsub(paste("(.*\\s*\\",locus_tag,"=)(.*)(\")", sep=""),"\\2", 
                      grep(locus_tag,entry, value=T )[1])
      try(z[z$id==i, "old_locus_tag"]<-
            gsub("(.*\\s*\\old_locus_tag=\")(.*)(\")","\\2", 
                 grep("old_locus_tag=",entry, value=T )[1]), silent=T)
      
      z[z$id==i,"locus_tag"]<-
        gsub("(.*?\")(.*)", "\\2",  pre_locus)
      z[z$id==i, "direction"]<-
        ifelse(!any(grepl("complement", entry_head)), "leading", "compliment")
      try(z[z$id==i, "note"]<-
            gsub("  ","",paste0(entry[grep("note=", entry): 
                                        (grep("codon_start",entry)-1)], collapse="  ")),
          silent=T)
      
    } else if(any(grepl("tRNA", entry_head))){
      z[z$id==i, "type"]<-"tRNA"
      pre_locus<-gsub(paste("(.*\\s*\\",locus_tag,"=)(.*)(\")", sep=""),"\\2", 
                      grep(locus_tag,entry, value=T )[1])
      try(z[z$id==i, "old_locus_tag"]<-
            gsub("(.*\\s*\\old_locus_tag=\")(.*)(\")","\\2", 
                 grep("old_locus_tag=",entry, value=T )[1]), silent=T)
      
      z[z$id==i,"locus_tag"]<-
        gsub("(.*?\")(.*)", "\\2",  pre_locus)
      z[z$id==i, "direction"]<-
        ifelse(!any(grepl("complement", entry_head)), "leading", "compliment")
      try(z[z$id==i, "product"]<-
            gsub("(.*\\s*/product=)(\")(.*)(\")","\\3", grep("product",entry, value=T )),
          silent=T)
      try(z[z$id==i, "inference"]<-
            gsub("(.*\\s*\\inference=)(*\")(.*)","\\3",
                 grep("inference",entry, value=T )), silent=T)
      try(z[z$id==i, "anticodon"]<-
            gsub("(.*\\s*/anticodon=)(.*)","\\2", 
                 grep("anticodon",entry, value=T )),silent=T)
      #####
      #   cant use the normal entry_head because occasionally other entries will 
      #   say "rRNA" in header
    } else if(any(grepl("rRNA", entry[c(1,entry_end)]))){
      z[z$id==i, "type"]<-"rRNA"
      pre_locus<-gsub(paste("(.*\\s*\\",locus_tag,"=)(.*)(\")", sep=""),"\\2", 
                      grep(locus_tag,entry, value=T )[1])
      try(z[z$id==i, "old_locus_tag"]<-
            gsub("(.*\\s*\\old_locus_tag=\")(.*)(\")","\\2", 
                 grep("old_locus_tag=",entry, value=T )[1]), silent=T)
      
      z[z$id==i,"locus_tag"]<-
        gsub("(.*?\")(.*)", "\\2",  pre_locus)
      z[z$id==i, "direction"]<-
        ifelse(!any(grepl("complement", entry_head)), "leading", "compliment")
      try(z[z$id==i, "note"]<-
            gsub("  ","",paste0(entry[grep("note=", entry): 
                                        (grep("codon_start",entry)-1)], collapse="  ")),
          silent=T)
      z[z$id==i, "product"]<-
        gsub("(.*\\s*/product=)(\")(.*)(\")","\\3", grep("product",entry, value=T ))
      #####
    } else if(any(grepl("misc_feature", entry_head))){
      z[z$id==i, "type"]<-"misc_feature"
      try(z[z$id==i, "old_locus_tag"]<-
            gsub("(.*\\s*\\old_locus_tag=\")(.*)(\")","\\2", 
                 grep("old_locus_tag=",entry, value=T )[1]), silent=T)
      
      z[z$id==i, "direction"]<-
        ifelse(!any(grepl("complement", entry_head)), "leading", "compliment")
      try(z[z$id==i, "note"]<-
        gsub("  ","",paste0(entry[grep("note=\"", entry): (length(entry)-1)], collapse="  ")),silent=T)
      #####
    } else if(any(grepl("misc_binding", entry_head))){
      z[z$id==i, "type"]<-"misc_binding"
      try(z[z$id==i, "old_locus_tag"]<-
            gsub("(.*\\s*\\old_locus_tag=\")(.*)(\")","\\2", 
                 grep("old_locus_tag=",entry, value=T )[1]), silent=T)
      
      z[z$id==i, "direction"]<-ifelse(!any(grepl("complement", entry_head)), 
                                      "leading", "compliment")
      try(z[z$id==i, "note"]<-
            gsub("  ","",
                 paste0(entry[grep("note=", entry):
                                (grep("codon_start",entry)-1)],collapse="  ")
            ), silent=T)
      try(z[z$id==i, "bound_moiety"]<-
            gsub("(.*\\s*/bound_moiety=\")(.*)(\")","\\2",
                 grep("bound_moiety",entry, value=T )), silent=T)
      #####
    } else if(any(grepl("assembly_gap", entry_head))){
      z[z$id==i, "type"]<-"assembly_gap"
      try(z[z$id==i, "old_locus_tag"]<-
            gsub("(.*\\s*\\old_locus_tag=\")(.*)(\")","\\2", 
                 grep("old_locus_tag=",entry, value=T )[1]), silent=T)
      
      try(z[z$id==i, "estimated_length"]<-
            gsub("(.*\\s*/estimated_length=)(.*)","\\2",
                 grep("estimated_length",entry, value=T )), silent=T)
      try(z[z$id==i, "gap_type"]<-
            gsub("(.*\\s*/gap_type=\")(.*)(\")","\\2",
                 grep("gap_type",entry, value=T )), silent=T)
      try(z[z$id==i, "linkage_evidence"]<-
            gsub("(.*\\s*/linkage_evidence=\")(.*)(\")","\\2",
                 grep("linkage_evidence",entry, value=T )), silent=T)
      ##### 
    } else if(any(grepl("repeat_region", entry_head))){
      z[z$id==i, "type"]<-"repeat_region"
      try(z[z$id==i, "old_locus_tag"]<-
            gsub("(.*\\s*\\old_locus_tag=\")(.*)(\")","\\2", 
                 grep("old_locus_tag=",entry, value=T )[1]), silent=T)
      
      z[z$id==i, "repeat_region"]<-
        gsub("  ","",paste0(grep("(\\d*\\.\\.\\d*)", entry[1], value=T), collapse="  "))
    } else if(any(grepl("misc_RNA", entry_head))){
      z[z$id==i, "type"]<-"misc_RNA"
      try(z[z$id==i, "product"]<-
            gsub("(.*\\s*/product=)(\")(.*)(\")","\\3", grep("product",entry, value=T )),
          silent=T)
      try(z[z$id==i, "old_locus_tag"]<-
            gsub("(.*\\s*\\old_locus_tag=\")(.*)(\")","\\2", 
                 grep("old_locus_tag=",entry, value=T )[1]), silent=T)
      try(z[z$id==i, "gene"]<-
            gsub("(.*\\s*\\gene=\")(.*)(\")","\\2", 
                 grep("gene=",entry, value=T )[1]), silent=T)
    } else if(any(grepl("ncRNA", entry_head))){
      z[z$id==i, "type"]<-"ncRNA"
      try(z[z$id==i, "product"]<-
            gsub("(.*\\s*/product=)(\")(.*)(\")","\\3", grep("product",entry, value=T )),
          silent=T)
      try(z[z$id==i, "old_locus_tag"]<-
            gsub("(.*\\s*\\old_locus_tag=\")(.*)(\")","\\2", 
                 grep("old_locus_tag=",entry, value=T )[1]), silent=T)
      try(z[z$id==i, "gene"]<-
            gsub("(.*\\s*\\gene=\")(.*)(\")","\\2", 
                 grep("gene=",entry, value=T )[1]), silent=T)
    } else if(any(grepl("tmRNA", entry_head))){
      z[z$id==i, "type"]<-"tmRNA"
      try(z[z$id==i, "product"]<-
            gsub("(.*\\s*/product=)(\")(.*)(\")","\\3", grep("product",entry, value=T )),
          silent=T)
      try(z[z$id==i, "old_locus_tag"]<-
            gsub("(.*\\s*\\old_locus_tag=\")(.*)(\")","\\2", 
                 grep("old_locus_tag=",entry, value=T )[1]), silent=T)
      try(z[z$id==i, "gene"]<-
            gsub("(.*\\s*\\gene=\")(.*)(\")","\\2", 
                 grep("gene=",entry, value=T )[1]), silent=T)
    } else if(any(grepl("STS", entry_head))){
      z[z$id==i, "type"]<-"STS"
      try(z[z$id==i, "db_xref"]<-
            gsub("(.*\\s*/db_xref=)(\")(.*)(\")","\\3", 
                 grep("db_xref",entry, value=T )), silent=T)
      try(z[z$id==i, "old_locus_tag"]<-
            gsub("(.*\\s*\\old_locus_tag=\")(.*)(\")","\\2", 
                 grep("old_locus_tag=",entry, value=T )[1]), silent=T)
      
      try(z[z$id==i, "standard_name"]<-
            gsub("(.*\\s*/standard_name=)(\")(.*)(\")","\\3", 
                 grep("standard_name",entry, value=T )), silent=T)
    } else if(any(grepl("source", entry_head))){
      z[z$id==i, "type"]<-"source"
      try(z[z$id==i, "old_locus_tag"]<-
            gsub("(.*\\s*\\old_locus_tag=\")(.*)(\")","\\2", 
                 grep("old_locus_tag=",entry, value=T )[1]), silent=T)
      
      z[z$id==i, "source"]<-
        gsub("  ","",paste0(entry[grep("organism=", entry): 
                                    entry_end], collapse="  "))
    }
  }
  z<-z[1:(nrow(z)-1),]    #  bye bye dummy row
  z
}
##  IMPORTANT!!!  THIS RELIES ON THE RANGES DF HAVING THE LAST ROW AS BUFFER
#       Why is this important?  if you wanna subset, you must include 
#        the last row, either by using tail() or rbind(subset, tail(gl1,1))
#         if you dont, you will get an error.  
# 

#^^^^^^^^^ 
#ranges<-extract_features(data=returns[[2]],ranges =  ranges,locus_tag = "locus_tag")

################################################################################
# this function extracts the nucleotide sequence from the .gb file, along with 
#upstream and downstream regions of interest
get_seqs<-function(gbdf, seq, upstream=500, downstream=500){
  idlist<-gbdf[gbdf$type !="source","id"]
  for( i in idlist){ #print(i)}
    #gbdf["id"==i, "seq"]<- 
    gbdf[gbdf$id==i,"dnaseq"]<-
      substr(seq, gbdf[gbdf$id==i, "loc_start"], gbdf[gbdf$id==i, "loc_end"])
    gbdf[gbdf$id==i,"dnaseq_upstream"]<-
      substr(seq, gbdf[gbdf$id==i, "loc_start"]-upstream, 
             gbdf[gbdf$id==i, "loc_start"]-1)
    gbdf[gbdf$id==i,"dnaseq_downstream"]<-
      substr(seq, gbdf[gbdf$id==i, "loc_end"], 
             gbdf[gbdf$id==i, "loc_end"]+downstream)
    if (i==idlist[length(idlist)]){
      next
    } else{
     # print(paste("id:",i))
      # 15 is the approx limit to the RBS
    gbdf[gbdf$id==i,"preceeding_intergenic_region"]<-
      ifelse(gbdf[gbdf$id==i, "direction"]=="leading",
        substr(seq, gbdf[gbdf$id==i-1, "loc_end"]+15,gbdf[gbdf$id==i, "loc_start"]-15),
        substr(seq, gbdf[gbdf$id==i, "loc_end"]+15,gbdf[gbdf$id==i+1, "loc_start"]-15))
    }
    #gbdf["id"==i, "upstrea"]
  }
  gbdf$region<-paste(gbdf$dnaseq_upstream, toupper(gbdf$dnaseq), gbdf$dnaseq_downstream, sep="")
  gbdf$search_region<-
    ifelse(gbdf$direction=="leading",
           paste(gbdf$dnaseq_upstream, toupper(substr(gbdf$dnaseq, 1, 50)), sep=""),
           paste(toupper(substr(gbdf$dnaseq, nchar(gbdf$dnaseq)-50, nchar(gbdf$dnaseq))), gbdf$dnaseq_downstream, sep=""))
  gbdf$prom_region<- 
    ifelse(gbdf$direction=="leading",
           paste(substr(gbdf$dnaseq_upstream,   1, (nchar(gbdf$dnaseq_upstream)-15)), sep=""),
           paste(substr(gbdf$dnaseq_downstream, 15, nchar(gbdf$dnaseq_downstream)), sep=""))

  gbdf
}
#^^^^^^^^^ 
#result<-get_seqs(gbdf = ranges, seq = returns[[3]], upstream = 500, downstream = 500)
##################  from seqinR, plus my edits
c2s<-function (chars = c("m", "e", "r", "g", "e", "d")) {
  return(paste(chars, collapse = ""))
}
# s2c<-function (string) {
#   if (is.character(string) && length(string) == 1) {
#     return(.Call("s2c", string, PACKAGE = "seqinr"))
#   }
#   else {
#     warning("Wrong argument type in s2c(), NA returned")
#     return(NA)
#   }
# }

write.fasta<-function (sequences, names, file.out, open = "w", nbchar = 60, 
                       as.string = FALSE) {
  outfile <- file(description = file.out, open = open)
  write.oneseq <- function(sequence, name, nbchar, as.string) {
    writeLines(paste(">", name, sep = ""), outfile)
    if (as.string) {
      #sequence2 <- s2c(sequence)
      seq<-strsplit(sequence, NULL)[[1]] #instead of s2c
    } else {
      seq<-sequence
    }
    #get dims
    l <- length(seq) #total chars
#    l <- nchar(sequence)   #  why length here?
    q <- floor(l/nbchar) # number of lines needed
    r <- l - nbchar * q # number of characters on last line
    if (q > 0) {
      lapply(1:q, function(x) {
        writeLines(c2s(seq[(nbchar * (x - 1) + 1):(nbchar * x)]), outfile)
      })
    }
    if (r > 0) {  #  write last line
      writeLines(c2s(seq[(nbchar * q + 1):l]), outfile)
    }
  }
  if (!is.list(sequences)) {
    write.oneseq(sequence = sequences, name = names, nbchar = nbchar, 
                 as.string = as.string)
  }
  else {
    n.seq <- length(sequences)
    lapply(1:n.seq, function(x){
       write.oneseq(sequence = as.character(sequences[[x]]),
                    name = names[x], nbchar = nbchar, as.string = as.string)
    })
  }
  close(outfile)
}
##################

################################################################################
#  Alrighty now, time to make this work for multiple sscaffolds. Test with a single
#  scaffold first, then move on to more ambitious pursuits of having this search the
#  folder for the scaf*.txt files
#  UPDATE this, as of 20151208, works with uams-1 genome with 2 scaffolds.

split_if_scaffolded(source.file)
print(paste("Found ",length(grep("scaf.*", dir())),
            " scaffold(s) in the current directory: (", getwd(), ")", sep=""))
print(grep("scaf.+", dir(), value=T))

scafList<-as.list(grep("scaf.*", dir(), value=T))
#for( i in grep("scaf.*", dir(), value=T)){
#resultsAll<-  lapply(scafList, function(i){
#fastaString<-""
resultsEachScaf<-  lapply(scafList, function(i){
  #i=scafList[2]
  #i=source.file
  i=unlist(i)
  outName<-gsub("(.*)(\\.txt)","\\1", i )
  returns<-load_gb(i)
  returns[[3]]<-clean_sequence(seq = returns[[3]])
  ###  give a fasta
#  seqs<-clean_sequence(seq = returns[[3]])
  #
  ranges<-get_ranges(returns[[2]])
  ranges<-clean_ranges_rna(ranges, returns[[2]])
  ranges<-extract_features(data=returns[[2]],ranges =  ranges,locus_tag = "locus_tag")
  result<-get_seqs(gbdf = ranges, seq = returns[[3]], upstream = upstream, downstream = downstream)
  result$scaffold<-gsub("(.*)(_scaf.*)","\\1", i)
  return(result)
  }
)
all_cols<-unique(unlist(lapply(resultsEachScaf, colnames))) #get all column names used
resultsEachScaf<-lapply(1:length(resultsEachScaf), function(x){
  if (length(all_cols[(!all_cols %in% colnames(resultsEachScaf[[x]]))])>0){
    new_cols<-all_cols[(!all_cols %in% colnames(resultsEachScaf[[x]]))]
    for(i in new_cols){
      resultsEachScaf[[x]][,i]<-NA
    }
  }
  return(resultsEachScaf[[x]])
})


#  extract just the sequences for scaffold fasta output
seqData<-  lapply(scafList, function(i){
  #i=scafList[2]
  #i=source.file
  i=unlist(i)
  outName<-gsub("(.*)(\\.txt)","\\1", i )
  returns<-load_gb(i)
  seqs<-clean_sequence(seq = returns[[3]])
  out<-list(i, seqs)
  return(out)
})
sequences<-list()
names<-list()
for (i in 1:length(seqData)){
  names<-c(names, gsub("(.*)(_scaf.*)","\\1",seqData[[i]][1]))
  sequences<-c(sequences,seqData[[i]][2])
}


input_name<-gsub("(.*)(\\.gb)","\\1", gsub("(.*)(\\/)(.*\\.gb)", "\\3", source.file))
#exclude<-c("misc_feature","misc_RNA","misc_binding","assembly_gap","source")  
exclude<-c("source","noSOFA", "assembly_gap")  

####  make a GFF3 file
if(length(resultsEachScaf)>1){
  finalDF<-  resultsEachScaf[[1]]
  attributes<-paste(
    "ID=",finalDF$type,":",finalDF$locus_tag,
    ";Name=", finalDF$gene, 
    ";Alias=", finalDF$old_locus_tag, 
    # ";Parent=",
    # ";Target=",
    # ";Gap=",
    # ";Derives_from=",
    # ";Note=",
    ";Dbxref=", finalDF$db_xref,
    # ";Ontology_term=", 
    # ";Is_circular=",
    sep='')
  attributesGTF<-paste(
    'gene_id "',finalDF$locus_tag,'";',
    'transcript_id "', finalDF$locus_tag, '";',
   sep='')
  #  give SOFA accession numbers for type
  gff3type<-lapply(finalDF$type, function(x){
    if (x =="CDS"){
      return("SO:0000316")
    } else if (x =="pseudo"){
      return("SO:0000336")
    } else if (x =="rRNA"){
      return("SO:0000252")
    } else if (x =="tRNA"){
      return("SO:0000253")
    } else if (x =="misc_feature"){ #genbank feature
      return("SO:0000253")
    } else if (x =="misc_RNA"){ #snoRNA
      return("SO:0000275")
    } else if (x =="misc_binding"){  #binding site
      return("SO:0000409")
   } else if (x =="assembly_gap"){ #  assembly component
      return("SO:0000143")
    }
  })
  gff3df<-data.frame("seqid"     = finalDF$scaffold,
                     "source"    ="Genbank via gbParse.R",
                     "type"      = finalDF$type,
                     "start"     = finalDF$loc_start,
                     "end"       = finalDF$loc_end,
                     "score"     = ".",
                     "strand"    = ifelse(finalDF$direction=="leading", "+","-"),
                     "phase"     = "0",
                     "attributes"= attributes)
  gtfdf<-data.frame("seqid"     = finalDF$scaffold,
                     "source"    ="Genbank via gbParse.R",
                     "type"      = finalDF$type,
                     "start"     = finalDF$loc_start,
                     "end"       = finalDF$loc_end,
                     "score"     = ".",
                     "strand"    = ifelse(finalDF$direction=="leading", "+","-"),
                     "frame"     = "0",
                     "attributes"= attributesGTF)
  mainHeaderA<-"##gff-version 3"
  scafHeaderA<-paste("##sequence-region",gff3df[1,"seqid"], 1,gff3df[1,"end"], sep=" ")
  gff3df<-gff3df[!(gff3df$type %in% exclude),] 
  gtfdf<-gtfdf[!(gtfdf$type %in% exclude),] 
  ####  hack up to add codons
  gtfdf <-gtfdf [!(gtfdf$type  %in% exclude),] 
  start_codons<-gtfdf[gtfdf$type=="CDS",]
  start_codons$type<-"start_codon"
  start_codons$end<-start_codons$start+2
  stop_codons<-gtfdf[gtfdf$type=="CDS",]
  stop_codons$type<-"stop_codon"
  stop_codons$start<-stop_codons$end-2
  exons<-gtfdf[gtfdf$type=="CDS",]
  exons$type<-"exon"
  new_cds<-gtfdf[gtfdf$type=="CDS",]
  new_cds$end <-  new_cds$end-2
  #  rip and stitch
  gtfdf<-gtfdf[gtfdf$type!="CDS",]
  gtfdf<-rbind(gtfdf,start_codons)
  gtfdf<-rbind(gtfdf,stop_codons)
  gtfdf<-rbind(gtfdf,exons)
  gtfdf<-rbind(gtfdf,new_cds)
  gtfdf<-gtfdf[order(gtfdf$start),]
  
  ######
  print(paste("writing gff file:", paste(working_dir, input_name, ".gff", sep="")))
  
  # writeLines(text=paste(c(mainHeaderA,scafHeaderA),sep = "\n"),
  #            con =  paste(working_dir, input_name, ".gff", sep=""))
  # write.table(x = gff3df, 
  #             append = T,
  #             col.names = F,row.names = F,
  #             file = paste(working_dir, input_name, ".gff", sep=""), 
  #             sep="\t", quote = F)
  gtfdf$strand<-ifelse( is.na(gtfdf$strand),".",as.character(gtfdf$strand))
  write.table(x = gtfdf, 
              append = F,
              col.names = F,row.names = F,
              file = paste(working_dir, input_name, ".gtf", sep=""), 
              sep="\t", quote = F)
######################################################################## 
  for( i in 2:length(resultsEachScaf)){
    interm<-resultsEachScaf[[i]]
    finalDF<-rbind(finalDF, interm)
    attributes<-paste(
      "ID=",interm$type,":",interm$locus_tag,
      ";Name=", interm$gene, 
      ";Alias=", interm$old_locus_tag, 
      # ";Parent=",
      # ";Target=", 
      # ";Gap=", 
      # ";Derives_from=", 
      # ";Note=", 
      ";Dbxref=", interm$db_xref,
      # ";Ontology_term=",
      # ";Is_circular=",
      sep='')
    attributesGTF<-paste(
      'gene_id "',interm$locus_tag,'";',
      'transcript_id "', interm$locus_tag, '";',
      sep='')
    gff3df<-data.frame("seqid"     = interm$scaffold,
                       "source"    ="Genbank via gbParse.R",
                       "type"      = unlist(interm),
                       "start"     = interm$loc_start,
                       "end"       = interm$loc_end,
                       "score"     = ".",
                       "strand"    = ifelse(interm$direction=="leading", "+","-"),
                       "phase"     = "0",
                       "attributes"= attributes)
    gtfdf<-data.frame("seqid"     = interm$scaffold,
                      "source"    ="Genbank via gbParse.R",
                      "type"      = interm$type,
                      "start"     = interm$loc_start,
                      "end"       = interm$loc_end,
                      "score"     = ".",
                      "strand"    = ifelse(interm$direction=="leading", "+","-"),
                      "frame"     = "0",
                      "attributes"= attributesGTF)
    #    mainHeaderA<-"##gff-version 3"
    scafHeaderA<-paste("##sequence-region",gff3df[1,"seqid"], 1,gff3df[1,"end"], sep=" ")
    gff3df<-gff3df[!(gff3df$type %in% exclude),] 
    gtfdf<-gtfdf[!(gtfdf$type %in% exclude),] 
    ####  hack up to add codons
    gtfdf <-gtfdf [!(gtfdf$type  %in% exclude),] 
    start_codons<-gtfdf[gtfdf$type=="CDS",]
    start_codons$type<-"start_codon"
    start_codons$end<-start_codons$start+2
    stop_codons<-gtfdf[gtfdf$type=="CDS",]
    stop_codons$type<-"stop_codon"
    stop_codons$start<-stop_codons$end-2
    exons<-gtfdf[gtfdf$type=="CDS",]
    exons$type<-"exon"
    new_cds<-gtfdf[gtfdf$type=="CDS",]
    new_cds$end <-  new_cds$end-2
    #  rip and stitch
    gtfdf<-gtfdf[gtfdf$type!="CDS",]
    gtfdf<-rbind(gtfdf,start_codons)
    gtfdf<-rbind(gtfdf,stop_codons)
    gtfdf<-rbind(gtfdf,exons)
    gtfdf<-rbind(gtfdf,new_cds)
    gtfdf<-gtfdf[order(gtfdf$start),]
    
    ######

    # cat(scafHeaderA,append = T,sep="\n",
    #     file =  paste(working_dir, input_name, ".gff", sep=""))
    # write.table(x = gff3df, 
    #             append = T,
    #             col.names = F,row.names = F,
    #             file = paste(working_dir, input_name, ".gff", sep=""), 
    #             sep="\t", quote = F)
    gtfdf$strand<-ifelse( is.na(gtfdf$strand),".",as.character(gtfdf$strand))
    write.table(x = gtfdf, 
                append = T,
                col.names = F,row.names = F,
                file = paste(working_dir, input_name, ".gtf", sep=""), 
                sep="\t", quote = F)
  }
} else{
  finalDF<-resultsEachScaf[[1]]
  attributes<-paste(
    "ID=",finalDF$type,":",finalDF$locus_tag,
    ";Name=", finalDF$gene, 
    ";Alias=", finalDF$old_locus_tag, 
    # ";Parent=",
    # ";Target=",
    # ";Gap=",
    # ";Derives_from=",
    # ";Note=",
    ";Dbxref=", finalDF$db_xref,
    # ";Ontology_term=", 
    # ";Is_circular=",
    sep='')
  attributesGTF<-paste(
    'gene_id "',finalDF$locus_tag,'";',
    'transcript_id "', finalDF$locus_tag, '";',
    sep='')
  gff3type<-lapply(finalDF$type, function(x){
    if (x =="CDS"){
      return("SO:0000316")
    } else if (x =="pseudo"){
      return("SO:0000336")
    } else if (x =="rRNA"){
      return("SO:0000252")
    } else if (x =="tRNA"){
      return("SO:0000253")
    } else if (x =="misc_feature"){ #genbank feature
      return("SO:0000253")
    } else if (x =="misc_RNA"){ #snoRNA
      return("SO:0000275")
    } else if (x =="misc_binding"){  #binding site
      return("SO:0000409")
    } else if (x =="assembly_gap"){ #  assembly component
      return("SO:0000143")
    } else {
      return("noSOFA")
    }
  })
  
  gff3df<-data.frame("seqid"     = finalDF$scaffold,
                     "source"    ="Genbank via gbParse.R",
                     "type"      = unlist(gff3type),
                     "start"     = finalDF$loc_start,
                     "end"       = finalDF$loc_end,
                     "score"     = ".",
                     "strand"    = ifelse(finalDF$direction=="leading", "+","-"),
                     "phase"     = "0",
                     "attributes"= attributes)
  gtfdf<-data.frame("seqid"     = finalDF$scaffold,
                    "source"    ="Genbank via gbParse.R",
                    "type"      = finalDF$type,
                    "start"     = finalDF$loc_start,
                    "end"       = finalDF$loc_end,
                    "score"     = ".",
                    "strand"    = ifelse(finalDF$direction=="leading", "+","-"),
                    "frame"     = "0",
                    "attributes"= attributesGTF, stringsAsFactors = F)
  mainHeaderA<-"##gff-version 3"
  scafHeaderA<-paste("##sequence-region",gff3df[1,"seqid"], 1,gff3df[1,"end"], sep=" ")
  gff3df<-gff3df[!(gff3df$type %in% exclude),]
####  hack up to add codons
  gtfdf <-gtfdf [!(gtfdf$type  %in% exclude),] 
  start_codons<-gtfdf[gtfdf$type=="CDS",]
  start_codons$type<-"start_codon"
  start_codons$end<-start_codons$start+2
  stop_codons<-gtfdf[gtfdf$type=="CDS",]
  stop_codons$type<-"stop_codon"
  stop_codons$start<-stop_codons$end-2
  exons<-gtfdf[gtfdf$type=="CDS",]
  exons$type<-"exon"
  new_cds<-gtfdf[gtfdf$type=="CDS",]
  new_cds$end <-  new_cds$end-2
  #  rip and stitch
  gtfdf<-gtfdf[gtfdf$type!="CDS",]
  gtfdf<-rbind(gtfdf,start_codons)
  gtfdf<-rbind(gtfdf,stop_codons)
  gtfdf<-rbind(gtfdf,exons)
  gtfdf<-rbind(gtfdf,new_cds)
  gtfdf<-gtfdf[order(gtfdf$start),]
  
######

  print(paste("writing gtf file:", paste(working_dir, input_name, ".gff", sep="")))
  
  # writeLines(text=paste(c(mainHeaderA,scafHeaderA),sep = "\n"),
  #            con =  paste(working_dir, input_name, ".gff", sep=""))
  # write.table(x = gff3df, 
  #             append = T,
  #             col.names = F,row.names = F,
  #             file = paste(working_dir, input_name, ".gff", sep=""), 
  #             sep="\t", quote = F)
  gtfdf$strand<-ifelse( is.na(gtfdf$strand),".",as.character(gtfdf$strand))
  write.table(x = gtfdf, 
              append = F,
              col.names = F,row.names = F,
              file = paste(working_dir, input_name, ".gtf", sep=""), 
              sep="\t", quote = F)
}

#  And they all lived happily ever after
for (fn in scafList){
  print(paste("removing temp file", fn))
if (file.exists(fn)) {file.remove(fn)}}
print(paste("writing csv file:",paste(working_dir, input_name, ".csv", sep="")))

write.csv(finalDF, paste(working_dir, input_name, ".csv", sep=""))


# 
# # writeLines(fastaString,
# #            con = paste(working_dir, input_name, ".fasta", sep=""))
print(paste("writing fasta file:",paste(working_dir, input_name, ".fasta", sep="")))
write.fasta(names = names,
            sequences =sequences,as.string = T,
              nbchar = 80,
            file.out =  paste(working_dir, input_name, ".fasta", sep=""))


##


