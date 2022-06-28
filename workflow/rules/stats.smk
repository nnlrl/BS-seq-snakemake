# # ==========================================================================================
# # Alignment stats extracted with samtools
#
# are included into multiqc report

rule idxstats_samtools:
    input:
        DIR_sorted+"{prefix}.bam",
    output:
        DIR_sorted+"{prefix}.idxstats.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        nice("samtools" ,["idxstats {input} > {output}"])


rule stat_samtools:
    input:
        DIR_sorted+"{prefix}.bam",
    output:
        DIR_sorted+"{prefix}.stats.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        nice("samtools", ['stats {input} > {output}'])


rule flagstat_samtools:
    input:
        DIR_sorted+"{prefix}.bam",
    output:
        DIR_sorted+"{prefix}.flagstat.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        nice("samtools" ,["flagstat {input} > {output}"])


# ==========================================================================================
# Create a csv file tabulating the lengths of the chromosomes in the reference genome:

rule tabulate_seqlengths:
    input:
        GENOMEFILE
    output:
        index       = GENOMEFILE+".fai",
        chrom_seqlengths  = os.path.join(DIR_mapped,ASSEMBLY+"_chromlengths.csv")
    message: fmt("Tabulating chromosome lengths in genome: {ASSEMBLY} for later reference.")
    conda:
        "../envs/samtools.yaml"
    shell:
        nice('samtools',
        ['faidx {input}',";",
        tool('cut'),"-f1,2","{output.index}","> {output.chrom_seqlengths}"])
