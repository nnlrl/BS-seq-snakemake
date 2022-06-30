## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      Create tabix files from methylDackel .methylKit files

      Arguments:
      --location location of input file
      --sample.id unique name of input sample
      --assembly assembly used to map the reads
      --treatment treatment of input sample
      --context context of methylation
      --mincov minimum coverage (default: 10)
      --dbdir name of the output folder
      --logFile file to print the logs to
      --help              - print this text

      Example:
      ./test.R --arg1=1 --arg2='output.txt' --arg3=TRUE \n\n")



  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")

argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))

names(argsL) <- argsDF$V1


## catch output and messages into log file
if(!is.null(argsL$logFile)) {
  out <- file(argsL$logFile, open = "at")
  sink(out,type = "output")
  sink(out, type = "message")
}



# Run Functions -----------------------------------------------------------


### Create tabix file from Methylation Calls

## load methylKit
# require("methylKit")

location  <- argsL$location
sample.id <- argsL$sample.id
assembly  <- argsL$assembly
treatment <- argsL$treatment
context   <- argsL$context
mincov    <- as.numeric(argsL$mincov)
dbdir     <- argsL$dbdir

message("read file <",location,"> into methylKit object")

## read file into methylKit object
methRaw = methylKit::methRead(location = location,
                    sample.id = sample.id,
                    assembly = assembly,
                    treatment = treatment,
                    mincov = mincov,
                    context = context,
                    dbtype = 'tabix',
                    dbdir = dbdir)

