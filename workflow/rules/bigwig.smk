# ==========================================================================================
# Export methylation from tabix to bigwig
#

def destrand(context):
    return config['general']['export-bigwig']['context'][context.lower()]['destrand']

rule export_tabix_bigwig:
    input:
        seqlengths = os.path.join(
            DIR_mapped, ASSEMBLY + "_chromlengths.csv"),
        filepath = os.path.join(
            DIR_methcall, "{tool}", "tabix_{context}/{prefix}_{context}.txt.bgz")
    output:
        bw = os.path.join(DIR_bigwig, "{prefix}.{context}_{tool}.bw")
    wildcard_constraints:
        context=".+(?<!destranded)"
    params:
        assembly = ASSEMBLY,
        destrand = lambda wc: destrand(wc.context)
    log:
        DIR_bigwig + "{prefix}.{context}.{tool}.export_tbx2bw.log"
    message: fmt("exporting methylation as bigwig for sample {wildcards.prefix}.")
    conda:
        "../envs/Rscript.yaml"
    shell:
        nice('Rscript', ["{DIR_scripts}/export_tbx2bw.R",
                         "--filepath={input.filepath}",
                         "--seqlengths_path={input.seqlengths}",
                         "--assembly={params.assembly}",
                         "--destrand={params.destrand}",
                         "--out_path={output}",
                         "--logFile={log}"], "{log}")


rule export_tabix_bigwig_destrand:
    input:
        seqlengths = os.path.join(
            DIR_mapped, ASSEMBLY + "_chromlengths.csv"),
        filepath = os.path.join(
            DIR_methcall, "{tool}", "tabix_{context}/{prefix}_{context}.txt.bgz")
    output:
        bw = os.path.join(DIR_bigwig, "{prefix}.{context}_destranded_{tool}.bw")
    params:
        assembly = ASSEMBLY,
        destrand = lambda wc: destrand(wc.context)
    log:
        DIR_bigwig + "{prefix}.{context}.{tool}.export_tbx2bw.log"
    message: fmt("exporting methylation as bigwig for sample {wildcards.prefix}.")
    conda:
        "../envs/Rscript.yaml"
    shell:
        nice('Rscript', ["{DIR_scripts}/export_tbx2bw.R",
                         "--filepath={input.filepath}",
                         "--seqlengths_path={input.seqlengths}",
                         "--assembly={params.assembly}",
                         "--destrand={params.destrand}",
                         "--out_path={output}",
                         "--logFile={log}"], "{log}")
