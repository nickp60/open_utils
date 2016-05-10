########################################################################
#                            NargParse.R
#                              version 0.2.3
#                             by Nick Waters 20160421
# minor version changes:
#  improved arg_bool to have defaults
test_args<-c("~/GitHub/R/")
arg_handle_help<-function(args_vector,  version="0.0", def_args=test_args, help_message2=help_message){
  # elipsis is to pass help_message down function chain
  if(is.null(help_message2)){ help_message2="help message here"}
  CONTINUE_bool=F
  y=" " # begin with space to separate collapsed args
  if(length(args_vector)==0) stop("No arguments given!  run with -h to see usage")
  for (i in 1:length(args_vector)){
    y<-paste(y,args_vector[i], sep=" ")
  }
  y=paste(y," ", sep="") #add spaces to collapsed args
  if (!is.vector(y)){
    stop("args not a vector!")
  }
  if(!is.character(help_message2)){
    stop("invalid help message")
  }
  if(grepl(" -H | -h | --help ",y)){
    cat(help_message2)
    stop("Stopping script. Rerun without help flag")
  }
  if(grepl(" -v | -V | --version ",y)){
    cat(paste("version:",version))
    stop("Stopping script. Rerun without version flag")
  }
  CONTINUE_bool=T
  return(CONTINUE_bool)
}
#############################  Positional Arg checks #################
arg_directory<-function(x, need_write=T, VERBOSE=T, makedir=F, recursive=F){
  if(is.na(x)) stop("man, I thought I was the only one to pass NA's to this function!" )
  if(file.access(x, mode = 0)==0 & dir.exists(x)){
    ifelse(VERBOSE,print("directory exists"),"")
    if(file.access(x, mode = 4)==0 & need_write==F){
      ifelse(VERBOSE,print("directory is readable"),"")
    } else if ((file.access(x, mode = 4)==0 & need_write==T) &
               (file.access(x, mode=2)!=0)){
      stop("directory is readable but this script requires writing permissions")
    } else if (file.access(x, mode=2)==0){
      ifelse(VERBOSE,print("directory is readable and writable"),"")
    }
    return_val<-x
  } else{
    if(makedir){
      hosting_dir<-dirname(x)
      if (file.access(hosting_dir, mode = 2)==0  & dir.exists(hosting_dir)){
        ifelse(VERBOSE, print("creating new directory"),'')
        ifelse(VERBOSE, x,'')
        dir.create(path = x, showWarnings = T, recursive = recursive)
        return_val<-x
      }
    } else {
      return_val<-NULL
      stop(paste(x, "is not a valid directory"))
    }
  }
  return(return_val)
}

arg_file<-function(x){
  if(!file.exists(x)){
    stop("this isnt a valid file")
  }
  return(x)
}
arg_string<-function(x, to_lower=F){
  if(dir.exists(x) | file.exists(x)){
    stop("this is the path to a file or directory; cannot be used!")
  } else if (is.character(x)){
  #TODO
  }
  if(to_lower){
    x=tolower(x)
  }
  return(x)
}
arg_integer<-function(x){
  try(x<-as.integer(x))
  return(x)
}
arg_numeric<-function(x, type=c( "numeric")){
  if(type=="integer"){
    try(x<-as.integer(x))
    return(x)
  } else if (type=="numeric"){
    try(x<-as.numeric(x))
    return(x)
  }
}
arg_bool<-function(x, default=F){
  if(is.numeric(x)) stop("cant convert numeric arguments to logical!")
  if(is.na(as.logical(x))){
    return(default)
  } else{
  return(as.logical(x))
  }
}
######################   Flagged Args ################################
#  this will be agnostic to things like '-h' vs '--help'.  Too much work for now.
#  this will be for pythonic flags with the dashes, but not exclusively, to
#    avoid the issue with negative numbers as arguments

#  The output wll be a named vector without the dashes 

#  x is your commandArgs vector
# this will check to make sure all the required flags are there, and parse all to named vector without the dashes
parse_flagged_args<-function(x, required=NULL, optional=NULL, help_message=NULL,
                             version="0.0",test_args=test_args,DEBUG=F){
  if(DEBUG){
    x=c("-i", "~/GitHub/R/","-j", "15","-o", "none")
    required=c("-i","-j" )
    optional=c("-o")
  }
  if(!arg_handle_help(args_vector = x, version =  version, help_message2 = help_message)){
    stop("Try again without version or help flags")
  } else{
    if(DEBUG) print("parsing args...")
    }
  options_used<-c(required, optional[optional%in%x])
  for( i in options_used){if(length(grep(paste0("^",i,"$"), x))>1) stop("Do you have repeated flags?")}
  args_output<-c()
  # check against numeric flags
  if(!all(is.na(suppressWarnings(sapply(c(optional, required), as.numeric))))){ stop("Cannot have numeric flags!")}
  # ensure at least one required
  if(length(required)<1){
    stop("must have at least one required argument, at this point anyway...")
  }
  for(i in 1:(length(options_used))){
    #check for missing requiements
    if(options_used[i] %in% required & !any(grepl(options_used[i], x))){
      stop(paste("Required arguments: ", required, ";\n missing ",options_used[i], sep="" ))
    }
    x_loc<-grep(paste(options_used[i],"$", sep=''), x)
    thing<-x[x_loc+1] # get the next thing right after the flag
    #
    names(thing)<-gsub('^-', '', options_used[i])
    args_output<-c(args_output,thing)
  }
  return(args_output)
}



