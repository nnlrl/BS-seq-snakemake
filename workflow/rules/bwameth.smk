# # # ==========================================================================================
# # # Mapping bisulfite reads with bwa-meth:
# #

# Notes:
# bwa meth paper https://arxiv.org/pdf/1401.1129.pdf
# https://www.biostars.org/p/93726/
# bwa meth and snp http://seqanswers.com/forums/showthread.php?t=40562

#
# # ==========================================================================================
# # Generate methyl-converted version of the reference genome:
#

rule bwameth_genome_preparation:
    input:
        ancient(GENOMEFILE)
    output:
        GENOMEFILE+".bwameth.c2t.sa",
        GENOMEFILE+".bwameth.c2t.amb",
        GENOMEFILE+".bwameth.c2t.ann",
        GENOMEFILE+".bwameth.c2t.pac",
        GENOMEFILE+".bwameth.c2t.bwt",
        GENOMEFILE+".bwameth.c2t"
    log:
        os.path.join(GENOMEPATH,'bwameth_genome_preparation.output_'+ASSEMBLY+'.log')
    message: "Converting {ASSEMBLY} Genome into Bisulfite analogue with bwa-meth"
    conda:
        "../envs/bwameth.yaml"
    shell:
        nice("bwameth", ["index {input}"],"{log}")


rule bwameth_touch_index:
    input:
        GENOMEFILE+".bwameth.c2t.sa",
        GENOMEFILE+".bwameth.c2t.amb",
        GENOMEFILE+".bwameth.c2t.ann",
        GENOMEFILE+".bwameth.c2t.pac",
        GENOMEFILE+".bwameth.c2t.bwt"
    output:
        GENOMEFILE+".bwameth.c2t_was.touched"
    message: "Update timestamp for {ASSEMBLY} Genome Index"
    shell:
        "sleep 60; touch {input};touch {output}"


# # ==========================================================================================
# # Align whole data sets
#


def bwameth_input(sample):
    files = list_files_TG(samplesheet(sample,'files'),sample, '')
    return(files)

rule bwameth_align_trimmed:
    input:
        rules.bwameth_touch_index.output,
        index = rules.bwameth_genome_preparation.output,
        files = lambda wc: bwameth_input(wc.sample)
    output:
        bam = DIR_mapped+"{sample}.bwameth.bam"
    log:
        DIR_mapped+"{sample}_bwameth_mapping.log"
    message: fmt("Mapping reads to genome using bwa-meth for sample {wildcards.sample}.")
    threads:
        config['execution']['rules']['bwameth_align_trimmed']['threads']
    conda:
        "../envs/bwameth.yaml"
    shell:
      nice("bwameth",["--reference {GENOMEFILE}","-t {threads}",
          "{input.files}","2> {log}","|",
          tool('samtools'), "view -bS - > {output.bam}"])


# ==========================================================================================
# Deduplication by samblaster
#
# https://github.com/GregoryFaust/samblaster
#
rule samblaster_markdup_sort:
    input:
        DIR_mapped+"{sample}.bwameth.bam"
    output:
        bam = DIR_sorted+"{sample}.bwameth.sorted.markdup.bam",
        index = DIR_sorted+"{sample}.bwameth.sorted.markdup.bam.bai"
    params:
        memory=config['execution']['rules']['samblaster_markdup_sort']['memory']
    log:
        DIR_sorted+"{sample}_markdups.log"
    message: fmt("Deduplicating reads with samblaster for sample {wildcards.sample}")
    threads:
        config['execution']['rules']['samblaster_markdup_sort']['threads']
    conda:
        "../envs/samtools.yaml"
    shell:
        nice("samtools",
        ["view -h {input}"," | ",
        tool("samblaster"),toolArgs("samblaster"),"2>> {log}","|",
        tool("samtools"),"sort",
         "-o {output.bam}", "-@ {threads}",
         "-m {params.memory}", "-l 9","2>> {log}",";",
         tool("samtools"),"index {output.bam}"],("{log}"))
