## Wrapper function to run a specific Rmd script
## which does the computation and generates an HTML report with code included
render2HTML <- function(reportFile,
                            outFile,
                            outDir,
                            finalReportDir,
                            report.params=NULL,
                            self.contained=TRUE,
                            logFile = NULL)
{

  # knitr 1.39 changed the default behavior.
  options(knitr.graphics.rel_path = FALSE)

  output_format = rmarkdown::all_output_formats(reportFile, 'UTF-8')

  if(!dir.exists(finalReportDir)) dir.create(finalReportDir, recursive = TRUE)

  if(is.null(report.params)) report.params <- list()

  ## render single report
  message("render single report")

  ## make independent intermediate dirs
  interDir <- paste0(outDir,"/",outFile,"_tmp")

  htmlwidgets::setWidgetIdSeed(1234)
  rmarkdown::render(
    input = reportFile,
    output_dir = outDir,
    intermediates_dir = interDir,
    output_file = outFile,
    knit_root_dir = outDir,
    output_format = rmarkdown::html_document(
      toc = TRUE,
      toc_float = TRUE,
      theme = 'lumen',
      number_sections = FALSE,
      code_folding = "hide",
      self_contained = self.contained,
      includes = list(in_header = "pigx_bsseq_logo.html"),
      bibliography= "reports.bib"
    ),
    params=report.params,
    quiet = FALSE,
    clean = TRUE
  )


  on.exit(unlink(interDir,recursive = TRUE),add = TRUE)


}

