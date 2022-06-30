# ========================================================================================
# Process bam files into methyl-called formats:

rule bam_methCall:
    input:
        bamfile         = os.path.join(DIR_sorted,"{prefix}.bam")
    output:
        tabixfile       = os.path.join(DIR_methcall,"methylKit","tabix_{context}","{prefix}_{context}.txt.bgz"),
        tabixfileindex  = os.path.join(DIR_methcall,"methylKit","tabix_{context}","{prefix}_{context}.txt.bgz.tbi")
    params:
        assembly        = ASSEMBLY,
        mincov          = int(config['general']['methylation-calling']['minimum-coverage']),
        minqual         = int(config['general']['methylation-calling']['minimum-quality']),
        context         = "{context}"
    log:
        os.path.join(DIR_methcall,"{prefix}_{context}_meth_calls.log")
    message: fmt("Extract methylation calls from bam file.")
    conda:
        "../envs/Rscript.yaml"
    threads:
        int(config["execution"]["rules"]["methCall"]["threads"])
    shell:
        nice('Rscript', ["{DIR_scripts}/methCall.R",
                         "--inBam={input.bamfile}",
                         "--assembly={params.assembly}",
                         "--cores={threads}",
                         "--mincov={params.mincov}",
                         "--minqual={params.minqual}",
                         "--context={params.context}",
                         "--tabix={output.tabixfile}",
                         "--logFile={log}"],"{log}")


# =========================================================================================
# Perform segmentation on the methylome:

rule methseg:
    ## paths inside input and output should be relative
    input:
        tabixfile   =  os.path.join(DIR_methcall,"{tool}","tabix_{context}","{prefix}_{context}.txt.bgz")
    output:
        bedfile     = os.path.join(DIR_seg,"{prefix}_{context}_{tool}.meth_segments.bed")
    params:
        methSegPng  = os.path.join(DIR_seg,"{prefix}_{context}_{tool}.meth_segments.png"),
        sample_id   = "{prefix}",
        assembly    = ASSEMBLY
    log:
        os.path.join(DIR_seg,"{prefix}_{context}_{tool}.meth_segments.log")
    message: fmt("Segmenting methylation profile for {wildcards.context} context of sample {wildcards.prefix}.")
    conda:
        "../envs/Rscript.yaml"
    shell:
        nice('Rscript', ["{DIR_scripts}/methSeg.R",
                         "--tabix={input.tabixfile}",
                         "--outBed={output.bedfile}",
                         "--png={params.methSegPng}",
                         "--sample.id={params.sample_id}",
                         "--assembly={params.assembly}",
                         "--logFile={log}"],"{log}")
