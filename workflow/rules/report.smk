rule final_report:
    input:
        methCall_tabixfile  = os.path.join(DIR_methcall,"{tool}","tabix_{context}","{prefix}_{context}.txt.bgz"),
        methSeg_bedfile     = os.path.join(DIR_seg,"{prefix}_{context}_{tool}.meth_segments.bed"),
        template            = os.path.join(DIR_templates,"index.Rmd"),
        bigwigFile          = os.path.join(DIR_bigwig, "{prefix}.{context}_{tool}.bw")
    output:
        report        = os.path.join(DIR_final,
                                     "sample_reports",
                                     "{prefix}_{context}_{tool}_{assembly}_final.html")
    params:
        ## absolute path to bamfiles
        sampleid                     = "{prefix}",
        source_dir                   = config['locations']['input-dir'],
        out_dir                      = OUTDIR,
        inBam                        = os.path.join(OUTDIR, DIR_sorted,"{prefix}.sorted.markdup.bam"),
        context                      = "{context}",
        assembly                     = ASSEMBLY,
        mincov                       = int(config['general']['methylation-calling']['minimum-coverage']),
        minqual                      = int(config['general']['methylation-calling']['minimum-quality']),
        TSS_plotlength               = int(config['general']['reports']['TSS_plotlength']),
        methSegPng                   = os.path.join(DIR_seg,"{prefix}_{context}_{tool}.meth_segments.png"),
        scripts_dir                  = DIR_scripts,
        refGene_bedfile              = REFGENES_BEDFILE,
        webfetch                     = WEBFETCH,
        # required for any report
        bibTexFile                   = BIBTEXPATH,
        prefix                       = "{prefix}_{context}_{tool}_{assembly}_{tool}",
        workdir                      = os.path.join(DIR_final,"sample_reports"),
    log:
        os.path.join(DIR_final,"sample_reports", "{prefix}_{context}_{tool}_{assembly}_final.log")
    message: fmt("Compiling final report for sample {wildcards.prefix}.")
    # run:
    #     generateReport(input, output, params, log, "")
    conda:
        "../envs/Rscript.yaml"
    shell:
        nice('Rscript', ["{DIR_scripts}/generate_report.R",
                           "--reportFile={input.template}",
                           "--outFile={output.report}",
                           "--workdir={params.workdir}",
                           "--bibTexFile={params.bibTexFile}",
                           "--prefix={params.prefix}",
                           "--report.params='{{"+
                           ",".join([
                               '"sampleid":"{params.sampleid}"',
                               '"assembly":"{params.assembly}"',
                               '"context":"{params.context}"',
                               '"bigwigFile":"{input.bigwigFile}"',
                               '"inBam":"{params.inBam}"',
                               '"methCall_tabixfile":"{input.methCall_tabixfile}"',
                               '"methSegBed":"{input.methSeg_bedfile}"',
                               '"methSegPng":"{params.methSegPng}"',


                               '"mincov":"{params.mincov}"',
                               '"minqual":"{params.minqual}"',
                               '"TSS_plotlength":"{params.TSS_plotlength}"',

                               '"source_dir":"{params.source_dir}"',
                               '"out_dir":"{params.out_dir}"',
                               '"scripts_dir":"{params.scripts_dir}"',

                               '"refGenes_bedfile":"{params.refGene_bedfile}"',
                               '"webfetch":"{params.webfetch}"'
                           ])+"}}'",
                           "--logFile={log}"], "{log}",
                           "echo '' ")
