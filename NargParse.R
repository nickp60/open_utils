########################################################################
#                            NargParse.R
#                              version 0.1
#                             by Nick Waters 20160403
#   Simple functions for commandline tool development in R to provide simple ways of setting 
#      acceptable sanatized inputs
test_args<-c("~/GitHub/R/")

arg_handle_help<-function(x, help_message="help message here", version="0.0", test_args=test_args){
  y=" " # begin with space to separate collapsed args
  for (i in 1:length(x)){
    y<-paste(y,x[i], sep=" ")
  }
  y=paste(y," ", sep="") #add spaces to collapsed args
  if (!is.vector(y)){
    stop("args not a vector!")
  }
  if(!is.character(help_message)){
    stop("invalid help message")
  }
  if(grepl(" -H | -h | --help ",y)){
    print(help_message)
    print("such as:")
    print(test_args)
    stop("Stopping script. Rerun without help flag")
  }
  if(grepl(" -v | -V | --version ",y)){
    print(paste("version:",version))
    stop("Stopping script. Rerun without version flag")
  }
}

arg_directory<-function(x, need_write=T, VERBOSE=T){
  if(file.access(x, mode = 0)==0 & dir.exists(x)){
    ifelse(VERBOSE,print("directory exists"),"")
    if(file.access(x, mode = 4)==0 & need_write==F){
      ifelse(VERBOSE,"directory is readable","")
    } else if ((file.access(x, mode = 4)==0 & need_write==T) &
               (file.access(x, mode=2)!=0)){
      stop("directory is readable but this script requires writing permissions")
    } else if (file.access(x, mode=2)==0){
      ifelse(VERBOSE,"directory is readable and writable","")
    }
  } else{
    stop(paste(x, "is not a valid directory"))
  }
  return(x)
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
arg_numeric<-function(x){
  try(x<-as.numeric(x))
  return(x)
}
arg_bool<-function(x){
  try(x<-as.logical(x))
  return(x)
}