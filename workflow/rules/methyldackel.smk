
# ==========================================================================================
# Extract methylation counts with methylDackel
#


def protocol(wc):
    return config["SAMPLES"][wc.sample]['Protocol'].upper()


def keepDups(protocol):
    keepDupsFlag = str(config["general"]["methylation-calling"]["keep-duplicates"]).lower()
    keepDups = ""
    if keepDupsFlag == 'auto':
        if (protocol == "RRBS"):
            keepDups = "--keepDupes"
    elif keepDupsFlag in {'true','yes'}:
        keepDups = "--keepDupes"

    return keepDups


## TODO: compress methylKit files
## TODO: make extraction of chg and chh files optional??
rule methyldackel_extract_methylKit:
    input:
        bamfile = DIR_sorted + "{sample}.bwameth.sorted.markdup.bam",
        genome = GENOMEFILE
    output:
        cpgCallFile = DIR_methcall + "methylDackel/" +"{sample}_methyldackel_CpG.methylKit",
        chgCallFile = DIR_methcall + "methylDackel/" + "{sample}_methyldackel_CHG.methylKit",
        chhCallFile = DIR_methcall + "methylDackel/" + "{sample}_methyldackel_CHH.methylKit"
    wildcard_constraints:
        sample=".+(?<!deduped)"
    params:
        prefix = DIR_methcall + "methylDackel/" + "{sample}_methyldackel",
        protocol = lambda wc: protocol(wc),
        keepDups = lambda wc: keepDups(protocol(wc)),
        minqual = int(config['general']
                      ['methylation-calling']['minimum-quality'])
    log:
        DIR_methcall + "{sample}.methyldackel_calls.log"
    message: fmt("Extract methylation calls from bam file using MethylDackel for sample {wildcards.sample} and protocol {params.protocol}")
    threads:
        config['execution']['rules']['methyldackel_extract']['threads']
    conda:
        "../envs/methyldackel.yaml"
    resources:
        mem_mb = config["execution"]["rules"]["methyldackel_extract"]["memory"]
    shell:
        nice("methyldackel",
             ["extract", "{input.genome}", "{input.bamfile}",
              "-o {params.prefix}", "-@ {threads}", "{params.keepDups}",
              "--methylKit", "--CHH", "--CHG", "-q {params.minqual}"],
             ("{log}"))


## TODO: compress methylKit files
rule methyldackel_extract_methylKit_deduped:
    input:
        bamfile = DIR_sorted + "{sample}.bwameth.sorted.markdup.bam",
        genome = GENOMEFILE
    output:
        cpgCallFile = DIR_methcall + "methylDackel/" + \
            "{sample}.deduped_methyldackel_CpG.methylKit",
        chgCallFile = DIR_methcall + "methylDackel/" + \
            "{sample}.deduped_methyldackel_CHG.methylKit",
        chhCallFile = DIR_methcall + "methylDackel/" + \
            "{sample}.deduped_methyldackel_CHH.methylKit"
    params:
        prefix = DIR_methcall + "methylDackel/" + "{sample}.deduped_methyldackel",
        protocol = lambda wc: protocol(wc),
        keepDups = lambda wc: keepDups(protocol(wc)),
        minqual = int(config['general']
                      ['methylation-calling']['minimum-quality'])
    log:
        DIR_methcall + "{sample}.deduped.methyldackel_calls.log"
    message: fmt("Extract methylation calls from bam file using MethylDackel for sample {wildcards.sample} and protocol {params.protocol}")
    threads:
        config['execution']['rules']['methyldackel_extract']['threads']
    conda:
        "../envs/methyldackel.yaml"
    resources:
        mem_mb = config["execution"]["rules"]["methyldackel_extract"]["memory"]
    shell:
        nice("methyldackel",
             ["extract", "{input.genome}", "{input.bamfile}",
              "-o {params.prefix}", "-@ {threads}", "{params.keepDups}",
              "--methylKit", "--CHH", "--CHG", "-q {params.minqual}"],
             ("{log}"))


# ==========================================================================================
# Extract methylation bias with methylDackel
#

rule methyldackel_mbias:
    input:
        bamfile = DIR_sorted + "{sample}.bwameth.sorted.markdup.bam",
        genome = GENOMEFILE
    output:
        txt = DIR_methcall + "methylDackel/" + "{sample}_mbias_methyldackel.txt",
        p1 = DIR_methcall + "methylDackel/"+ "{sample}_mbias_OB.svg",
        p2 = DIR_methcall + "methylDackel/" + "{sample}_mbias_OT.svg"
    params:
        prefix = DIR_methcall + "methylDackel/" + "{sample}_mbias",
        protocol = lambda wc: protocol(wc),
        keepDups = keepDups("{params.protocol}"),
        minqual = int(config['general']
                      ['methylation-calling']['minimum-quality'])
    log:
        DIR_methcall + "{sample}.methyldackel_mbias.log"
    message: fmt("Calculate methylation bias using MethylDackel for sample {wildcards.sample}.")
    threads:
        config['execution']['rules']['methyldackel_extract']['threads']
    conda:
        "../envs/methyldackel.yaml"
    shell:
        nice("methyldackel",
             ["mbias", "{input.genome}", "{input.bamfile}",
              "{params.prefix}", "{params.keepDups}",
              "-@ {threads}", "--txt",
              "-q {params.minqual}", "> {output.txt}",
              "2> {log}"])


# ==========================================================================================
# Extract methylation bias with methylDackel
#

rule methyldackel_cytosine_report:
    input:
        bamfile = DIR_sorted + "{sample}.bwameth.sorted.markdup.bam",
        genome = GENOMEFILE
    output:
        DIR_methcall + "methylDackel/" + "{sample}_methyldackel.cytosine_report.txt"
    params:
        prefix = DIR_methcall + "methylDackel/" + "{sample}_methyldackel",
        protocol = lambda wc: protocol(wc),
        keepDups = keepDups("{params.protocol}"),
        minqual = int(config['general']
                      ['methylation-calling']['minimum-quality'])
    log:
        DIR_methcall + "{sample}.methyldackel_cytosine_report.log"
    message: fmt("Extract cytosine report from bam file using MethylDackel for sample {{sample}} and context {params.context}")
    threads:
        config['execution']['rules']['methyldackel_extract']['threads']
    conda:
        "../envs/methyldackel.yaml"
    resources:
        mem_mb = config["execution"]["rules"]["methyldackel_extract"]["memory"]
    shell:
        nice("methyldackel",
             ["extract", "{input.genome}", "{input.bamfile}",
              "-o {params.prefix}", "{params.keepDups}", "-@ {threads}",
              "--cytosine_report", "--CHH", "--CHG", "-q {params.minqual}"],
             ("{log}"))


# ==========================================================================================
# Convert to tabix files with methylKit
#

rule tabix_methyldackelfile:
    input:
        DIR_methcall + "methylDackel/" + "{prefix}_methyldackel_{context}.methylKit",
    output:
        DIR_methcall + "methylDackel/" + "tabix_{context}/{prefix}_{context}.txt.bgz",
        DIR_methcall + "methylDackel/" + "tabix_{context}/{prefix}_{context}.txt.bgz.tbi"
    params:
        sampleid = "{prefix}_{context}",
        assembly = ASSEMBLY,
        treatment = lambda wc: samplesheet(
            wc.prefix.replace(".deduped", ""), 'Treatment'),
        context = "{context}",
        dbdir = DIR_methcall + "methylDackel/" + "tabix_{context}/",
        mincov = int(config['general']['methylation-calling']
                     ['minimum-coverage'])
    log:
        DIR_methcall + "methylDackel/" + "tabix_{context}/{prefix}_{context}.makeTabix.log"
    message: fmt("Create Tabix file from MethylDackel file for sample {wildcards.prefix} and context {params.context}")
    conda:
        "../envs/Rscript.yaml"
    shell:
        nice('Rscript', ["{DIR_scripts}/makeTabix.R",
                         "--location={input}",
                         "--sample.id={params.sampleid}",
                         "--assembly={params.assembly}",
                         "--treatment={params.treatment}",
                         "--context={params.context}",
                         "--mincov={params.mincov}",
                         "--dbdir={params.dbdir}",
                         "--logFile={log}"], "{log}")
