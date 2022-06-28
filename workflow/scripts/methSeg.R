## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      Segment methylation profile using methylKit

      Arguments:
      --tabix     name of the input tabix file containting the methylRaw object
      --outBed    name of output BED file containing Segments
      --png       name of file to save diagnostic plots to
      --sample.id sample name
      --assembly  genome assembly
      --logFile   file to print the logs to
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

## Segmentation

## load methylKit
suppressPackageStartupMessages(library("methylKit"))

input     <- argsL$tabix
output    <- argsL$outBed
pngFile   <- argsL$png
sample.id <- argsL$sample.id
assembly  <- argsL$assembly


message("Reading tabix file.")
## read input tabix to methylRawDB
methRawDB <- methRead(location=input,
                      sample.id =sample.id,
                      assembly = assembly,
                      dbtype ="tabix")

# ## catch a possible error and touch empty files
# ## to trigger successful run
# err <- tryCatch(
#   expr = {
    ## try to run the code
    # png(filename = pngFile,
    #     units = "in",width = 8,
    #     height = 4.5,res=300)
    pdf(pngFile, 20, 20)

    message("Performing segmentation...")
    ### Segmentation of methylation profile
    res.gr = methSeg(methRawDB,
                     diagnostic.plot=TRUE)

    dev.off()

    ### Export

    message("Exporting segmentation...")
    ## export segments to bed file
    methSeg2bed(segments = res.gr,
                trackLine = paste0("track name='meth segments ' ",
                                   "description='meth segments of ",
                                   methRawDB@sample.id,
                                   " mapped to ",
                                   methRawDB@assembly,
                                   "' itemRgb=On"),
                colramp=colorRamp(c("gray","green", "darkgreen")),
                filename = output)
#   },
#   error = function(x) {
#     ## if it fails still generate empty output
#     file.create(output)
#     message(paste("error occured!!",x))
#   }
# )



