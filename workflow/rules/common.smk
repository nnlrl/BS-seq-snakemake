# ==============================================================================
#
#                                       HELPER FUNCTIONS
#
# put here helper functions that are used within more than one rule
# ==============================================================================

import os
import time
from glob import glob
from datetime import datetime

# --------------------------------------
# general purpose
# --------------------------------------

def fmt(message):
    """Format the MESSAGE string."""
    return "----------  " + message + "  ----------"

def tool(name):
    return config['tools'][name]['executable']

def toolArgs(name):
    if 'args' in config['tools'][name]:
        return config['tools'][name]['args']
    else:
        return ""

# sample sheet accessor functions
def samplesheet(name, item=None):
    """Access the SAMPLES dict from config."""
    if item:
        return config["SAMPLES"][name][item]
    else:
        config["SAMPLES"][name]

# print datetime next to a message
def log_time(text):
    _time = '$(date +"[%Y-%m-%d %T]")'
    # _time = time.strftime("%Y-%m-%d %T",time.localtime(time.time()))
    log = '{} {}'.format(_time, text)
    return(fmt(log))

# Generate a command line string that can be passed to snakemake's
# "shell".  The string is prefixed with an invocation of "nice".
def nice(cmd, args, log=None, fallback=None):
    executable = tool(cmd)
    line = ["nice -" + str(config['execution']['nice']),
            executable] + [toolArgs(cmd)] + args
    if log:
        line.insert(0, "echo {} > {};".format(log_time("Starting Now"),log))
        line.append(">> {} 2>&1".format(log))
    if fallback:
        line.append(" || {} ".format(fallback))
    if log:
        line.append("; echo {} >> {};".format(log_time("Done"),log))
    return " ".join(line)


# abandone current execution with a helpful error message:
def bail(msg):
    """Print the error message to stderr and exit."""
    print(msg, file=sys.stderr)
    exit(1)

def TrueOrFalse(value):
    value = repr(value)
    if value:
        answer = value.lower() in ["true","yes"]
    else:
        answer = False
    return(answer)


# check for common input/configuration errors:
def validate_config(config):
    # Check that all locations exist
    for loc in config['locations']:
        if ( (not loc == 'output-dir') and (not (os.path.isdir(config['locations'][loc]) or os.path.isfile(config['locations'][loc])))):
            bail("ERROR: The following necessary directory/file does not exist: {} ({})".format(
                config['locations'][loc], loc))

    # Check that all of the requested differential methylation
    # treatment values are found in the sample sheet.
    treatments = set([config["SAMPLES"][sample]["Treatment"]
                      for sample in config["SAMPLES"]])
    if 'DManalyses' in config:
        if config['DManalyses']:
            for analysis in config['DManalyses']:
                    for group in config['DManalyses'][analysis]['treatment_sample_groups'].split(",") + config['DManalyses'][analysis]['control_sample_groups'].split(","):
                        group = group.strip() #remove any leading/trailing whitespaces in the sample group names
                        if not any(treat == group for treat in treatments):
                            bail("ERROR: Invalid treatment group '{}' in analysis '{}'".format(
                            group, analysis))
        else:
            bail("ERROR: The config file contains empty 'DManalyses' section, please consider removing or commenting out this section.\n")


    if 'treatment-groups' in config['general']['differential-methylation']:
        bail("ERROR: The specification of treatment groups and differential analysis has changed.\n"+
        "Please retrieve the new default settings layout with 'pigx-bsseq --init settings'.\n")

    # Check for a any Assembly string
    if not config['general']['assembly']:
            bail("ERROR: Please set a genome assembly string in the settings file at general::assembly.")

    # Check for a any Assembly string
    if not (config['general']['use_bwameth'] or config['general']['use_bismark']):
            bail("ERROR: Please enable one or both bisulfite aligners at general::use_bwameth/use_bismark.")


    # Check if we have permission to write to the reference-genome directory ourselves
    # if not, then check if the ref genome has already been converted
    genome_dir = os.path.dirname(config['locations']['genome-fasta'])
    if (not os.access(genome_dir, os.W_OK) and
            not os.path.isdir(os.path.join(genome_dir, 'Bisulfite_Genome'))):
        bail("ERROR: reference genome has not been bisulfite-converted, and PiGx does not have permission to write to that directory. Please either (a) provide Bisulfite_Genome conversion directory yourself, or (b) enable write permission in {} so that PiGx can do so on its own.".format(
            genome_dir))

    # Check for a genome fasta file
    fasta = glob(os.path.join(genome_dir, '*.fasta'))
    fa    = glob(os.path.join(genome_dir, '*.fa'))
    if not len(fasta) + len(fa) == 1 :
        bail("ERROR: Missing (or ambiguous) reference genome: The number of files ending in either '.fasta' or '.fa' in the following genome directory does not equal one: {}".format(genome_dir))


# --------------------------------------
# sample related
# --------------------------------------

def isRRBS(args):
    sample = args[0]
    if config['SAMPLES'][sample]['Protocol'] == "RRBS":
        return "--rrbs"
    elif config['SAMPLES'][sample]['Protocol'] == "WGBS":
        return ""
    else:
        print("Protocol must be RRBS/WGBS")
        exit(1)

def dedupe_tag(protocol):
    if protocol.upper() == "WGBS":
        return ".deduped"
    elif protocol.upper() == "RRBS":
        return ""
    else:
        raise Exception("=== ERROR: unexpected protocol ===")


def get_fastq_name(full_name):
    # single end
    find_se_inx = full_name.find('_se_bt2')
    # paired-end
    find_pe_inx = full_name.find('_1_val_1_bt2')

    if(find_se_inx >= 0):
        output = full_name[:find_se_inx]
    elif(find_pe_inx >= 0):
        output = full_name[:find_pe_inx]
    else:
        bail("Unable to infer sample fastq name; cannot find trimming string in filename. \nHave the files been trimmed for adapter sequence and quality?")

    return(output)


# --------------------------------------
# context related
# --------------------------------------
def destrand(context):
    return config['general']['export-bigwig']['context'][context.lower()]['destrand']


# --------------------------------------
# generate dynamic output files per rule
# --------------------------------------

def files_for_sample(proc):
    return [expand(proc(config['SAMPLES'][sample]['files'],
                        config['SAMPLES'][sample]['SampleID'],
                        config['SAMPLES'][sample]['Protocol']))
            for sample in config['SAMPLES']]


def list_files_rawQC(files, sampleID, protocol):
    PATH = DIR_rawqc
    if len(files) == 1:
        return [PATH+sampleID+"_fastqc.zip"]  # ---- single end
    elif len(files) == 2:
        # ---- paired end
        return [PATH+sampleID+"_1_fastqc.zip", PATH+sampleID+"_2_fastqc.zip"]
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_TG(files, sampleID, protocol):
    PATH = DIR_trimmed
    if len(files) == 1:
        return [PATH+sampleID+"_trimmed.fq.gz"]  # ---- single end
    elif len(files) == 2:
        # ---- paired end
        return [PATH+sampleID+"_1_val_1.fq.gz", PATH+sampleID+"_2_val_2.fq.gz"]
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_posttrim_QC(files, sampleID, protocol):
    PATH = DIR_posttrim_QC
    if len(files) == 1:
        return [PATH+sampleID+"_trimmed_fastqc.zip"]  # ---- single end
    elif len(files) == 2:
        # ---- paired end
        return [PATH+sampleID+"_1_val_1_fastqc.zip", PATH+sampleID+"_2_val_2_fastqc.zip"]
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_bismark(files, sampleID, protocol):
    PATH = DIR_mapped
    if len(files) == 1:
        return [PATH+sampleID+"_trimmed_bismark_bt2_SE_report.txt",
                PATH+sampleID+"_trimmed_bismark_bt2.bam"]  # ---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_val_1_bismark_bt2_PE_report.txt",
                PATH+sampleID+"_1_val_1_bismark_bt2_pe.bam"]  # ---- paired end
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_bwameth(files, sampleID, protocol):
    PATH = DIR_mapped
    if len(files) == 1:
        return [PATH+sampleID+".bwameth.bam"]  # ---- single end
    elif len(files) == 2:
        return [PATH+sampleID+".bwameth.bam"]  # ---- paired end
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_dedupe(files, sampleID, protocol):
    PATH = DIR_sorted
    if len(files) == 1:
        # ---- single end
        return [PATH+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + ".bam"]
    elif len(files) == 2:
        # ---- paired end
        return [PATH+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + ".bam"]
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_markdup(files, sampleID, protocol):
    PATH = DIR_sorted
    if len(files) == 1:
        return [PATH+sampleID+".bwameth.sorted.markdup.bam"]  # ---- single end
    elif len(files) == 2:
        return [PATH+sampleID+".bwameth.sorted.markdup.bam"]  # ---- paired end
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_bwamethMappingStats(files, sampleID, protocol):
    PATH = DIR_sorted
    if len(files) == 1:
        return [PATH+sampleID+".bwameth.sorted.markdup.idxstats.txt",
                PATH+sampleID+".bwameth.sorted.markdup.stats.txt",
                PATH+sampleID+".bwameth.sorted.markdup.flagstat.txt"]  # ---- single end
    elif len(files) == 2:
        return [PATH+sampleID+".bwameth.sorted.markdup.idxstats.txt",
                PATH+sampleID+".bwameth.sorted.markdup.stats.txt",
                PATH+sampleID+".bwameth.sorted.markdup.flagstat.txt"]  # ---- paired end
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")



def bam_processing(files, sampleID, protocol):
    PATH = DIR_methcall+ "methylKit/"
    if len(files) == 1:
        # ---- single end
        return [PATH+"tabix_cpg/"+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + "_cpg.txt.bgz",
                PATH+"tabix_cpg/"+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + "_cpg.txt.bgz.tbi",
                PATH+"tabix_chg/"+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + "_chg.txt.bgz",
                PATH+"tabix_chg/"+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + "_chg.txt.bgz.tbi",
                PATH+"tabix_chh/"+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + "_chh.txt.bgz",
                PATH+"tabix_chh/"+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + "_chh.txt.bgz.tbi"
                ]
    elif len(files) == 2:
        # ---- paired end
        return [PATH+"tabix_cpg/"+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + "_cpg.txt.bgz",
                PATH+"tabix_cpg/"+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + "_cpg.txt.bgz.tbi",
                PATH+"tabix_chg/"+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + "_chg.txt.bgz",
                PATH+"tabix_chg/"+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + "_chg.txt.bgz.tbi",
                PATH+"tabix_chh/"+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + "_chh.txt.bgz",
                PATH+"tabix_chh/"+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + "_chh.txt.bgz.tbi"
                ]
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_methyldackel_extract(files, sampleID, protocol):
    PATH = DIR_methcall + "methylDackel/"
    return [PATH+sampleID+dedupe_tag(protocol)+"_methyldackel_CpG.methylKit",
            PATH+sampleID+dedupe_tag(protocol)+"_methyldackel_CHG.methylKit",
            PATH+sampleID+dedupe_tag(protocol)+"_methyldackel_CHH.methylKit",
            # PATH+sampleID+"_mbias_methyldackel.txt",
            # PATH+sampleID+"_mbias_OB.svg",
            # PATH+sampleID+"_mbias_OT.svg",
            # PATH+sampleID+"_methyldackel.cytosine_report.txt"
            ]


def list_files_maketabix_methyldackel(files, sampleID, protocol):
    PATH = DIR_methcall + "methylDackel/"
    return [
            PATH+"tabix_CpG/"+sampleID+dedupe_tag(protocol)+"_CpG.txt.bgz",
            PATH+"tabix_CpG/"+sampleID+dedupe_tag(protocol)+"_CpG.txt.bgz.tbi",
            PATH+"tabix_CHG/"+sampleID+dedupe_tag(protocol)+"_CHG.txt.bgz",
            PATH+"tabix_CHG/"+sampleID+dedupe_tag(protocol)+"_CHG.txt.bgz.tbi",
            PATH+"tabix_CHH/"+sampleID+dedupe_tag(protocol)+"_CHH.txt.bgz",
            PATH+"tabix_CHH/"+sampleID+dedupe_tag(protocol)+"_CHH.txt.bgz.tbi"
            ]


def bigwig_exporting_bismark(files, sampleID, protocol):
    PATH = DIR_bigwig
    res = []
    for context in ["cpg", "chg", "chh"]:
        DESTRAND = "_destranded" if destrand(context) else ""
        if len(files) == 1:
            # ---- single end
            res.append(PATH+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + "." + context + DESTRAND + "_methylKit"+ ".bw")
            # return PATH+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + ".cpg" + DESTRAND + "_methylKit"+ ".bw"
        elif len(files) == 2:
            # ---- paired end
            res.append(PATH+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + "." + context + DESTRAND +  "_methylKit" + ".bw")
            # return [PATH+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + ".cpg" + DESTRAND +  "_methylKit" + ".bw"]
        else:
            raise Exception(
                "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")
    return res

def bigwig_exporting_bwameth(files, sampleID, protocol):
    PATH = DIR_bigwig

    res = []
    for context in ["CpG", "CHH", "CHG"]:
        DESTRAND = "_destranded" if destrand(context.lower()) else ""
        if len(files) == 1:
            # ---- single end
            res.append(PATH+sampleID+dedupe_tag(protocol) + "." + context + DESTRAND + "_methylDackel" + ".bw")
            # return [PATH+sampleID+dedupe_tag(protocol) + "." + context + DESTRAND + "_methylDackel" + ".bw"]
        elif len(files) == 2:
            # ---- paired end
            res.append(PATH+sampleID+dedupe_tag(protocol) + "." + context + DESTRAND + "_methylDackel" + ".bw")
            # return [PATH+sampleID+dedupe_tag(protocol) + "." + context + DESTRAND + "_methylDackel" + ".bw"]
        else:
            raise Exception(
                "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")
    return res


def methSeg_bismark(files, sampleID, protocol):
    PATH = DIR_seg

    res = []
    for context in ["cpg", "chg", "chh"]:
        if len(files) == 1:
            # ---- single end
            res.append(PATH+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + "_" + context + "_methylKit" +".meth_segments.bed")
            # return PATH+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + "_cpg" + "_methylKit" +".meth_segments.bed"
        elif len(files) == 2:
            # ---- paired end
            res.append(PATH+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + "_" + context + "_methylKit" +".meth_segments.bed")
            # return [PATH+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + "_cpg" + "_methylKit" +".meth_segments.bed"]
        else:
            raise Exception(
                "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")
    return res


def methSeg_bwameth(files, sampleID, protocol):
    PATH = DIR_seg

    res = []
    for context in ["CpG", "CHG", "CHH"]:
        if len(files) == 1:
            # ---- single end
            return res
        elif len(files) == 2:
            # ---- paired end
            res.append(PATH+sampleID+dedupe_tag(protocol) + "_" + context + "_methylDackel" +".meth_segments.bed")
            # return [PATH+sampleID+dedupe_tag(protocol) + "_CpG" + "_methylDackel" +".meth_segments.bed"]
        else:
            raise Exception(
                "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")
    return res


def list_final_reports_bwameth(files, sampleID, protocol):
    res = []
    for context in ["CpG", "CHG", "CHH"]:
        SUFFIX = "{}_methylDackel_{}".format(context, ASSEMBLY)
        PATH = os.path.join(DIR_final,"sample_reports/" )
        if len(files) == 1:
            # ---- single end
            # return [PATH+sampleID+dedupe_tag(protocol) + ".CpG" + DESTRAND + "_methylDackel" + ".bw"]
            res.append(PATH+sampleID+dedupe_tag(protocol) + "_"+ SUFFIX+ "_final.html")
            # return [PATH+sampleID+dedupe_tag(protocol) + "_"+ SUFFIX+ "_final.html"]
        elif len(files) == 2:
            # ---- paired end
            res.append(PATH+sampleID+dedupe_tag(protocol) + "_"+ SUFFIX+ "_final.html")
            # return [PATH+sampleID+dedupe_tag(protocol) + "_"+ SUFFIX+ "_final.html"]
        else:
            raise Exception(
                "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")
    return res


def list_final_reports_bismark(files, sampleID, protocol):
    res = []
    for context in ["cpg", "chg", "chh"]:
        SUFFIX = "{}_methylKit_{}".format(context, ASSEMBLY)
        PATH = os.path.join(DIR_final,"sample_reports/" )
        if len(files) == 1:
            # ---- single end
            res.append(PATH+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + "_"+ SUFFIX+ "_final.html")
            # return PATH+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + "_"+ SUFFIX+ "_final.html"
        elif len(files) == 2:
            # ---- paired end
            res.append(PATH+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + "_"+SUFFIX+"_final.html")
            # return [PATH+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + "_"+SUFFIX+"_final.html"]
        else:
            raise Exception(
                "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")
    return res

# --------------------------------------
# diffmeth related
# --------------------------------------
def get_sampleids_from_treatment(treatment):
    """Get SampleIDs from treatment string."""
    sample_ids = list(config["SAMPLES"].keys())
    sample_treatments = [samplesheet(s,"Treatment") for s in sample_ids]

    sampleids_list = [sample_ids[i]
                 for i, x in enumerate(sample_treatments) if x == treatment]

    return(sampleids_list)


def get_sampleids_from_analysis(analysis):
    """Get SampleIDs for each Analysis group."""
    sampleids_list = []
    for group in config['DManalyses'][analysis]:
        for treatment in config['DManalyses'][analysis][group].split(","):
            sampleids_list += get_sampleids_from_treatment(treatment)

    return(sampleids_list)


def makeDiffMethPath(path, suffix, treatment):
    return path + str(treatment).replace('vs', '_') + dedupe_tag(config["SAMPLES"][get_sampleids_from_treatment(treatment[0])[0]]['Protocol']) + '_' + suffix


def diffmeth_input_function(treatments):
    treatments = treatments.replace(".deduped", "")
    sampleids = get_sampleids_from_treatment(treatments)

    inputfiles = []
    for sampleid in sampleids:
        fqname = config["SAMPLES"][sampleid]['fastq_name']
        protocol = config["SAMPLES"][sampleid]['Protocol']
        if len(fqname) == 1:
            inputfile = [os.path.join(
                DIR_methcall, sampleid+"_se_bt2.sorted" + dedupe_tag(protocol) + "_methylRaw.RDS")]
        elif len(fqname) == 2:
            inputfile = [os.path.join(
                DIR_methcall, sampleid+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + "_methylRaw.RDS")]
        inputfiles.append(inputfile)

    inputfiles = list(sum(inputfiles, []))
    return(inputfiles)


def files_for_treatment(proc):
    if "DManalyses" in config.keys():
        treatment_groups = config['DManalyses']
        if treatment_groups:
            return [expand(proc(comparison)) for comparison in treatment_groups if comparison]
        else:
            return []


def list_files_unite_bismark(treatment):
    PATH = DIR_diffmeth + treatment + "/"
    res = []
    context_list = config["DManalyses"][treatment]["context"].split(",")

    for context in context_list:
        DESTRAND = "_destranded" if destrand(context.lower()) else ""
        res.append(PATH + "methylBase_"+treatment+"_"+context.lower()+DESTRAND+"_methylKit.txt.bgz" )
    return res


def list_files_unite_bwameth(treatment):
    PATH = DIR_diffmeth + treatment + "/"
    res = []
    context_list = config["DManalyses"][treatment]["context"].split(",")

    for context in context_list:
        DESTRAND = "_destranded" if destrand(context.lower()) else ""
        res.append(PATH + "methylBase_"+treatment+"_"+context+DESTRAND+"_methylDackel.txt.bgz")
        # return [ PATH + "methylBase_"+treatment+"_"+context+DESTRAND+"_methylDackel.txt.bgz" ]
    return res


def list_files_diffmeth_bismark(treatment):
    PATH = DIR_diffmeth + treatment + "/"
    res = []
    context_list = config["DManalyses"][treatment]["context"].split(",")

    for context in context_list:
        DESTRAND = "_destranded" if destrand(context.lower()) else ""
        res.append(PATH + "methylDiff_"+treatment+"_"+context.lower()+DESTRAND+"_methylKit_full.txt.bgz")
        res.append(PATH + "methylDiff_"+treatment+"_"+context.lower()+DESTRAND+"_methylKit_results.tsv")
    return res
        # return [
        #         PATH + "methylDiff_"+treatment+"_cpg"+DESTRAND+"_methylKit_full.txt.bgz",
        #         PATH + "methylDiff_"+treatment+"_cpg"+DESTRAND+"_methylKit_results.tsv"
        #         ]


def list_files_diffmeth_bwameth(treatment):
    PATH = DIR_diffmeth + treatment + "/"
    res = []
    context_list = config["DManalyses"][treatment]["context"].split(",")

    for context in context_list:
        DESTRAND = "_destranded" if destrand(context.lower()) else ""
        res.append(PATH + "methylDiff_"+treatment+"_"+context+DESTRAND+"_methylDackel_full.txt.bgz")
        res.append(PATH + "methylDiff_"+treatment+"_"+context+DESTRAND+"_methylDackel_results.tsv")
    return res
        # return [
        #         PATH + "methylDiff_"+treatment+"_CpG"+DESTRAND+"_methylDackel_full.txt.bgz",
        #         PATH + "methylDiff_"+treatment+"_CpG"+DESTRAND+"_methylDackel_results.tsv"
        #         ]


def list_files_diffmeth_report_bwameth(treatment):
    PATH = DIR_final
    res = []
    context_list = config["DManalyses"][treatment]["context"].split(",")

    for context in context_list:
        DESTRAND = "_destranded" if destrand(context.lower()) else ""
        res.append(PATH + treatment + "/" + treatment + "_"+context+DESTRAND+"_methylDackel"+ ".diffmeth-report.html")
    return res
        # return [ PATH + treatment + "/" + treatment + "_CpG"+DESTRAND+"_methylDackel"+ ".diffmeth-report.html"]


def list_files_diffmeth_report_bismark(treatment):
    PATH = DIR_final
    res = []
    context_list = config["DManalyses"][treatment]["context"].split(",")

    for context in context_list:
        DESTRAND = "_destranded" if destrand(context.lower()) else ""
        res.append(PATH + treatment + "/" + treatment + "_"+context.lower()+DESTRAND+"_methylKit" +".diffmeth-report.html")
    return res
    # return [ PATH + treatment + "/" + treatment + "_cpg"+DESTRAND+"_methylKit" +".diffmeth-report.html"]
