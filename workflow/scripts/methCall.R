#methCall.R - takes a bismark bam file and exports methylKit tabix file



## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      Call methylation using methylKit

      Arguments:
      --inBam location of input bam file
      --assembly assembly used to map the reads
      --cores number of cores to use for methCall
      --mincov minimum coverage (default: 10)
      --minqual minimum base quality (default: 20)
      --context cytosine context to extract
      --tabix name of the tabix output file
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
out <- file(argsL$logFile, open = "at")
sink(out,type = "output")
sink(out, type = "message")



# Run Functions -----------------------------------------------------------

### Methylation Calling

## load methylKit
cores <- argsL$cores

suppressPackageStartupMessages(expr = {
library("methylKit")
data.table::setDTthreads(cores)
})


input     <- argsL$inBam
assembly  <- argsL$assembly
mincov    <- as.numeric(argsL$mincov)
minqual   <- as.numeric(argsL$minqual)
context   <- argsL$context
tabixfile   <- argsL$tabix

### Extract Methylation Calls

## extract the sample id from sample file
sample_id <- gsub(".bam","",basename(input))

## define the location to save intermediate file
save_folder <- dirname(tabixfile)

## format context string
contextStr <- switch(tolower(context),
                     "cpg" = "CpG",
                     "chg" = "CHG",
                     "chh" = "CHH"
                     )

message("Reading bam file into methylKit object")
## read bam file into methylKit object
methRawDB = processBismarkAln(location = input,
                            sample.id = sample_id,
                            assembly = assembly,
                            mincov = mincov,
                            minqual = minqual,
                            read.context = contextStr,
                            save.context = NULL,
                            save.folder = save_folder,
                            save.db=TRUE)

message("Tabix saved to: \n\t",getDBPath(methRawDB))

