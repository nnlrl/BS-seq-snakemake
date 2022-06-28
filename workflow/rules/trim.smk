# =========================================================================================
# Trim the reads for adapter-ends and quality

rule trim_reads_se:
    input:
       file = PATHIN + "{sample}.fq.gz"
    output:
       DIR_trimmed + "{sample}_trimmed.fq.gz" #---- this ALWAYS outputs .fq.qz format.
    params:
       extra      = config['tools']['trim_galore']['args'],
       outdir     = "--output_dir " + DIR_trimmed,
       phred      = "--phred33",
       gz         = "--gzip",
       isRRBS     = isRRBS
    log:
       DIR_trimmed + "{sample}.trimgalore.log"
    message: fmt("Trimming raw single-end read data from sample {wildcards.sample}")
    conda: "../envs/trim_galore.yaml"
    threads:
        config["execution"]["rules"]["trim_galore"]["threads"]
    shell:
       nice('trim_galore', ["-j {threads}", "{params}", "{input.file}"], "{log}")


rule trim_reads_pe:
    input:
        files = [ PATHIN+"{sample}_1.fq.gz",
                  PATHIN+"{sample}_2.fq.gz"]
    output:
        DIR_trimmed + "{sample}_1_val_1.fq.gz", #---- this ALWAYS outputs .fq.qz format.
        DIR_trimmed + "{sample}_2_val_2.fq.gz",
    params:
        extra          = config['tools']['trim_galore']['args'],
        outdir         = "--output_dir "+DIR_trimmed,
        phred          = "--phred33",
        gz             = "--gzip",
        paired         = "--paired",
        isRRBS         = isRRBS
    log:
        DIR_trimmed + "{sample}.trimgalore.log"
    message:
        fmt("Trimming raw paired-end read data from sample {wildcards.sample}")
    conda: "../envs/trim_galore.yaml"
    threads:
        config["execution"]["rules"]["trim_galore"]["threads"]
    shell:
        nice('trim_galore', ["-j {threads}", "{params}", "{input.files}"], "{log}")
