rule bismark_genome_preparation:
    input:
        file = GENOMEFILE,
        path = ancient(GENOMEPATH)
    output:
        GENOMEPATH+"Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
        GENOMEPATH+"Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"
    params:
        # pathToBowtie = "--path_to_bowtie "+ os.path.dirname(tool('bowtie2')),
        bismark_genome_preparation_args = config['tools']['bismark_genome_preparation']['args'],
        useBowtie2                      = "--bowtie2 ",
        verbose                         = "--verbose ",
        genomic_composition             = "--genomic_composition"
    log:
        os.path.join(GENOMEPATH,'bismark_genome_preparation_'+ASSEMBLY+'.log')
    message: fmt("Converting {ASSEMBLY} Genome into Bisulfite analogue")
    conda: "../envs/bismark.yaml"
    shell:
        nice('bismark_genome_preparation', ["{params}", "{input.path}"], "{log}")


bismark_cores = int(config['tools']['bismark']['cores'])


rule bismark_align_and_map_se:
    input:
        refconvert_CT = GENOMEPATH + "Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
	    refconvert_GA = GENOMEPATH + "Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa",
        fqfile = DIR_trimmed + "{sample}_trimmed.fq.gz",
        qc     = DIR_posttrim_QC + "{sample}_trimmed_fastqc.zip"
    output:
        DIR_mapped + "{sample}_trimmed_bismark_bt2.bam",
        DIR_mapped + "{sample}_trimmed_bismark_bt2_SE_report.txt"
    params:
        genomeFolder = "--genome_folder " + GENOMEPATH,
        outdir = "--output_dir  " + DIR_mapped,
        nucCov = "--nucleotide_coverage",
        useBowtie2   = "--bowtie2 ",
        tempdir      = "--temp_dir " + DIR_mapped,
        # cores = "--multicore " + bismark_cores
    log:
        DIR_mapped + "{sample}_bismark_se_mapping.log"
    message: fmt("Mapping single-end reads to genome {ASSEMBLY}")
    conda: "../envs/bismark.yaml"
    threads:
        bismark_cores
    shell:
        nice('bismark', ["--multicore {threads}", "{params}", "{input.fqfile}"], "{log}")


rule bismark_align_and_map_pe:
    input:
        refconvert_CT = GENOMEPATH + "Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
	    refconvert_GA = GENOMEPATH + "Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa",
        fin1 = DIR_trimmed + "{sample}_1_val_1.fq.gz",
        fin2 = DIR_trimmed + "{sample}_2_val_2.fq.gz",
        qc   = [ DIR_posttrim_QC + "{sample}_1_val_1_fastqc.zip",
                 DIR_posttrim_QC + "{sample}_2_val_2_fastqc.zip"]
    output:
        DIR_mapped + "{sample}_1_val_1_bismark_bt2_pe.bam",
        DIR_mapped + "{sample}_1_val_1_bismark_bt2_PE_report.txt"
    params:
        genomeFolder = "--genome_folder " + GENOMEPATH,
        outdir       = "--output_dir  " + DIR_mapped,
        nucCov       = "--nucleotide_coverage",
        useBowtie2   = "--bowtie2 ",
        tempdir      = "--temp_dir " + DIR_mapped,
        # cores        = "--multicore " + bismark_cores
    log:
        DIR_mapped+"{sample}_bismark_pe_mapping.log"
    message: fmt("Mapping paired-end reads to genome {ASSEMBLY}.")
    conda: "../envs/bismark.yaml"
    threads:
        bismark_cores
    shell:
        nice('bismark', ["--multicore {threads}", "{params}", "-1 {input.fin1}", "-2 {input.fin2}"], "{log}")


# =========================================================================================
# Deduplicate aligned reads from the bam file:

rule deduplication_se:
    input:
        DIR_mapped + "{sample}_trimmed_bismark_bt2.bam"
    output:
        DIR_sorted + "{sample}_se_bt2.sorted.deduped.bam"
    params:
        memory=config['execution']['rules']['samblaster_markdup_sort']['memory'],
        tmpdir=DIR_sorted+"{sample}/"
    log:
        DIR_sorted + "{sample}_deduplication.log"
    message: fmt("Deduplicating single-end aligned reads from {input}")
    conda: "../envs/samtools.yaml"
    threads:
        config['execution']['rules']['samblaster_markdup_sort']['threads']
    shell:
        nice("samtools",
        ["view -h {input}"," | ",
        "samblaster","-r", toolArgs("samblaster"),"2> {log}","|",
        "samtools","sort","-T={params.tmpdir}",
        "-o {output}", "-@ {threads}",
        "-m {params.memory}", "-l 9","2> {log}",";",
        "samtools","index {output}"],("{log}"))


#-----------------------
rule deduplication_pe:
    input:
        DIR_mapped + "{sample}_1_val_1_bismark_bt2_pe.bam"
    output:
        DIR_sorted + "{sample}_1_val_1_bt2.sorted.deduped.bam"
    params:
        memory=config['execution']['rules']['samblaster_markdup_sort']['memory'],
        tmpdir=DIR_sorted+"{sample}/"
    log:
        DIR_sorted+"{sample}_deduplication.log"
    message: fmt("Deduplicating paired-end aligned reads from {input}")
    conda: "../envs/samtools.yaml"
    threads:
        config['execution']['rules']['samblaster_markdup_sort']['threads']
    shell:
        nice("samtools",
        ["view -h {input}"," | ",
        "samblaster","-r",toolArgs("samblaster"),"2> {log}","|",
        "samtools","sort","-T={params.tmpdir}",
        "-o {output}", "-@ {threads}",
        "-m {params.memory}", "-l 9","2> {log}",";",
        "samtools","index {output}"],("{log}"))


# =========================================================================================
# Sort the bam file by position (and carry out mate-flagging in paired-end case):

rule sortbam_se:
    input:
        DIR_mapped+"{sample}_trimmed_bismark_bt2.bam"
    output:
        DIR_sorted+"{sample}_se_bt2.sorted.bam"
    message: fmt("Sorting bam file {input}")
    conda: "../envs/samtools.yaml"
    shell:
        nice('samtools', ["sort", "{input}", "-o {output}"])

#-----------------------
rule sortbam_pe:
    input:
        DIR_mapped+"{sample}_1_val_1_bismark_bt2_pe.bam"
    output:
        DIR_sorted+"{sample}_1_val_1_bt2.sorted.bam"
    message: fmt("Sorting bam file {input}")
    conda: "../envs/samtools.yaml"
    shell:
        nice('samtools', ["sort -n ", " {input} ", " | ", tool('samtools'), " fixmate -m  - - ", " | ", tool('samtools'), " sort -o {output} " ])
