# =========================================================================================
# Perform quality control on raw data

rule fastqc_raw: #----only need one: covers BOTH pe and se cases.
    input:
        PATHIN+"{sample}.fq.gz"
    output:
        DIR_rawqc+"{sample}_fastqc.zip"
    params:
        fastqc_args = config['tools']['fastqc']['args'],
        outdir      = "--outdir "+ DIR_rawqc     # usually pass params as strings instead of wildcards.
    log:
        DIR_rawqc+"{sample}_fastqc.log"
    message: fmt("Quality checking raw read data from sample {wildcards.sample}")
    conda: "../envs/fastqc.yaml"
    threads:
        config['execution']['rules']['fastqc']['threads']
    shell:
        nice('fastqc', ["-t {threads}", "{params}", "{input}"], "{log}")


# =========================================================================================
# Carry out post-trimming quality control

rule fastqc_after_trimming_se:
    input:
        DIR_trimmed+"{sample}_trimmed.fq.gz",
    output:
    	DIR_posttrim_QC+"{sample}_trimmed_fastqc.zip"
    params:
        fastqc_args = config['tools']['fastqc']['args'],
        outdir = "--outdir "+DIR_posttrim_QC
    log:
   	    DIR_posttrim_QC+"{sample}_trimmed_fastqc.log"
    message: fmt("Quality checking trimmmed single-end data from sample {wildcards.sample}")
    conda: "../envs/fastqc.yaml"
    threads:
        config['execution']['rules']['fastqc']['threads']
    shell:
        nice('fastqc', ["-t {threads}", "{params}", "{input}"], "{log}")


rule fastqc_after_trimming_pe:
    input:
        DIR_trimmed+"{sample}_1_val_1.fq.gz",
        DIR_trimmed+"{sample}_2_val_2.fq.gz"
    output:
    	DIR_posttrim_QC+"{sample}_1_val_1_fastqc.zip",
    	DIR_posttrim_QC+"{sample}_2_val_2_fastqc.zip"
    params:
        fastqc_args = config['tools']['fastqc']['args'],
        outdir = "--outdir "+DIR_posttrim_QC
    log:
   	    DIR_posttrim_QC+"{sample}_trimmed_fastqc.log"
    message: fmt("Quality checking trimmmed paired-end data from sample {wildcards.sample}")
    conda: "../envs/fastqc.yaml"
    threads:
        config['execution']['rules']['fastqc']['threads']
    shell:
        nice('fastqc', ["-t {threads}", "{params}", "{input}"], "{log}")


# ---------------------------------------------------------------------------- #
from itertools import chain

def multiqc_files(branch):
    files = []
    if config["general"]["trimming"]:
        files += files_for_sample(list_files_rawQC)
    files += files_for_sample(list_files_TG)
    files += files_for_sample(list_files_posttrim_QC)
    if branch == "bismark":
       files += files_for_sample(list_files_bismark)
    elif branch == "bwameth":
        files += files_for_sample(list_files_bwameth)
        files += files_for_sample(list_files_bwamethMappingStats)
    files = list(chain.from_iterable(files))
    return(files)


rule multiqc:
    input:
        files = lambda wc: multiqc_files(wc.branch)
    output:
        os.path.join(DIR_final,"multiqc", "{branch}_multiqc_report.html")
    params:
        trim_galore = DIR_trimmed,
        raw_qc = lambda wc: DIR_rawqc if config["general"]["trimming"] else "",
        posttrim_qc = DIR_posttrim_QC
    log:
        os.path.join(DIR_final,"multiqc", "{branch}_multiqc.log")
    message: "Generating multi-sample QC report for {wildcards.branch} branch."
    conda: "../envs/fastqc.yaml"
    shell:
      nice("multiqc",["-f","-n {output}","{input.files}","{params}"],"{log}")
