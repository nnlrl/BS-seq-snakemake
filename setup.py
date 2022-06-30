#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Author    : nnlrl
# Created   : 2022/05/27
# Title     : bsseq.py
# Objective :

description = """
BSseq Pipeline.

It is a data processing pipeline for raw fastq read data of
bisulfite experiments.  It produces methylation and coverage
information and can be used to produce information on differential
methylation and segmentation.
"""

import argparse
from os import path
import os, sys, json, csv, yaml, errno, shutil
from snakemake.utils import update_config
import shutil, filecmp
from glob import glob
import subprocess


WORKDIR = os.path.dirname(os.path.abspath(__file__))


def formatter(prog):
    return argparse.RawTextHelpFormatter(prog, max_help_position=80)


parser = argparse.ArgumentParser(description=description, formatter_class=formatter)

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument(
    "--init",
    dest="init",
    choices=["settings", "sample-sheet", "both"],
    const="both",
    nargs="?",
    help="Generate a template SETTINGS file, a SAMPLE-SHEET.  Leave empty for both.",
)

group.add_argument(
    "-s",
    "--start",
    dest="start",
    action="store_true",
    help="Start the pipeline using config/settings.yaml and config/sample_sheet.csv",
)

group.add_argument(
    "--clean",
    dest="clean",
    action="store_true",
    help="Clean up previous output files.",
)

parser.add_argument(
    "-c",
    "--configfile",
    dest="configfile",
    default="./config/config.json",
    help="""\
The config file used for calling the underlying snakemake process.  By
default the file 'config.json' is dynamically created from the sample
sheet and the settings file.
""",
)

parser.add_argument(
    "--target",
    dest="target",
    action="append",
    help="""\
Stop when the named target is completed instead of running the whole
pipeline.  The default target is "final-report". """,
)

parser.add_argument(
    "-n",
    "--dry-run",
    dest="dry_run",
    action="store_true",
    help="""\
Only show what work would be performed.  Do not actually run the
pipeline.""",
)

parser.add_argument(
    "--conda-create-envs-only",
    dest="conda_create_envs_only",
    action="store_true",
    help="""\
If specified, only creates the job-specific conda environments then exits. """,
)

parser.add_argument(
    "--graph",
    dest="graph",
    help="""\
Output a graph in PDF format showing the relations between rules of
this pipeline.  You must specify a graph file name such as
"graph.pdf".""",
)

parser.add_argument(
    "--force",
    dest="force",
    action="store_true",
    help="""\
Force the execution of rules, even though the outputs are considered
fresh.""",
)

parser.add_argument(
    "--reason",
    dest="reason",
    action="store_true",
    help="""\
Print the reason why a rule is executed.""",
)

parser.add_argument(
    "--unlock",
    dest="unlock",
    action="store_true",
    help="""\
Recover after a snakemake crash.""",
)

parser.add_argument(
    "--verbose",
    dest="verbose",
    action="store_true",
    help="""\
Print supplementary info on job execution.""",
)

args = parser.parse_args()


# Generate config file
def bail(msg):
    """Print the error message to stderr and exit."""
    print(msg, file=sys.stderr)
    exit(1)


def splitext_fqgz(string):
    def fq_suffix(filename):
        return any(filename.endswith(ext) for ext in [".fq", ".fastq", ".fasta"])

    def is_zipped(filename):
        return any(filename.endswith(ext) for ext in [".gz", ".bz2"])

    if is_zipped(string):
        string, zipext = os.path.splitext(string)
    else:
        zipext = ""
    if fq_suffix(string):
        base, ext = os.path.splitext(string)
        return (base, ext + zipext)
    else:
        bail("Input files are not fastq files!")


def get_filenames(mylist):
    return list(map(lambda x: splitext_fqgz(x)[0], mylist))


def parse_samples(lines):
    """
    Parse csv table with information about samples, eg:
    Read1,Read2,SampleID,ReadType,Treatment
    sampleB.pe1.fq.gz,sampleB.pe2.fq.gz,sampleB,WGBS,B,,
    pe1.single.fq.gz,,sampleB1,WGBS,B,,
    It returns a dictionary required for the config file.
    """
    sreader = csv.reader(lines, delimiter=",")
    all_rows = [row for row in sreader if row]

    header = list(map(lambda x: x.strip(), all_rows[0]))
    rows = all_rows[1:]
    minimal_header = ["Read1", "Read2", "SampleID", "Protocol", "Treatment"]

    if header[:5] != minimal_header:
        raise Exception(
            "First columns of the input table have to be "
            + ",".join(minimal_header)
            + "."
        )

    sample_ids = [x[2] for x in rows]
    if len(set(sample_ids)) != len(sample_ids):
        raise Exception("Column 'SampleID' has non-unique values.")

    # Create a dictionary with all params, keys are samples ids
    outputdict = {}
    for row in rows:
        if len(row) != 5:
            bail(
                "Invalid row format in Samplesheet. Each row should have five columns."
            )
        row = list(map(lambda x: x.strip(), row))
        files = list(filter(None, row[0:2]))
        if not files:
            raise Exception(
                "Each sample has to have an entry in at least one of the columns 'Read1' or 'Read2'."
            )

        sampleid_dict = {}
        for idx in range(len(header[2:])):
            try:
                sampleid_dict[header[2:][idx]] = row[2:][idx]
            except IndexError:
                raise Exception(
                    "Number of columns in row "
                    + idx
                    + " doesn't match number of elements in header."
                )

        sampleid_dict["files"] = files
        sampleid_dict["fastq_name"] = get_filenames(files)
        outputdict[row[2]] = sampleid_dict
    return {"SAMPLES": outputdict}


def generate_config(configfile, sample_sheet, settingsfile, dirs=None):
    """Generate a new configuration file CONFIGFILE using SAMPLE_SHEET and
    SETTINGSFILE as inputs.  Use the locations in DIRS to find default
    settings.
    """
    # Load defaults
    defaults = path.join(WORKDIR, "config", "settings.yaml")
    if not path.exists(defaults):
        bail(
            "Could not find default settings. You need to determine the file path of settings.yaml."
        )

    settings = yaml.safe_load(open(defaults, "r"))

    # Load user overrides
    if settingsfile:
        if not os.path.isfile(settingsfile):
            bail("ERROR: settings.yaml file not found")
        update_config(settings, yaml.safe_load(open(settingsfile, "r")))

    settings["execution"]["target"] = args.target

    # Record the location of the sample sheet.
    if not os.path.isfile(path.abspath(sample_sheet)):
        bail(
            "ERROR: Failed to find sample sheet provided. No such file: " + sample_sheet
        )
    settings["locations"]["sample-sheet"] = path.abspath(sample_sheet)

    # Load parameters specific to samples
    with open(sample_sheet, "r") as f:
        lines = f.read().splitlines()
    sample_params = parse_samples(lines)
    settings.update(sample_params)

    # add exec dir
    settings['locations'].update({"workdir": WORKDIR})

    # Resolve relative paths in the locations section
    # root = path.dirname(sample_sheet)
    # here = os.getenv('srcdir') if os.getenv('srcdir') else os.getcwd()

    for key in settings["locations"]:
        if settings["locations"][key]:
            settings["locations"][key] = path.normpath(
                path.join(WORKDIR, settings["locations"][key])
            )

    # Write the config file
    with open(configfile, "w") as outfile:
        dumps = json.dumps(
            settings,
            indent=4,
            sort_keys=True,
            separators=(",", ": "),
            ensure_ascii=True,
        )
        outfile.write(dumps)


# Create symbolic links to the inputs and reference genome

# Create links within the output folder that point directly to the
# reference genome, as well as to each sample input file so that it's
# clear where the source data came from.

# N.B. Any previously existing links will be kept in place, and no
# warning will be issued if this is the case.


def makelink(src, target):
    if not path.isfile(src):
        bail("Refusing to link non-existent file %s" % src)
    elif not path.isdir(path.dirname(target)):
        bail(
            "%s or subdirectory does not exist for linking %s"
            % config["locations"]["output-dir"],
            target,
        )
    else:
        try:
            os.symlink(src, target)
        except FileExistsError:
            pass


def maybe_copy(source, target_dir):
    target = path.join(target_dir, path.basename(source))
    if not path.exists(target) or not filecmp.cmp(source, target):
        shutil.copy(source, target_dir)


def prepare_links():
    os.makedirs(
        path.join(config["locations"]["output-dir"], "work"), exist_ok=True
    )

    os.makedirs(
        path.join(config["locations"]["output-dir"], "work/input"), exist_ok=True
    )

    # copy documentation file to the output work directory.
    contents = path.join(WORKDIR, "etc", "CONTENTS.txt")
    shutil.copyfile(
        contents, path.join(config["locations"]["output-dir"], "work/CONTENTS.txt")
    )

    # Link the reference genome
    try:
        os.symlink(
            path.dirname(config["locations"]["genome-fasta"]),
            path.join(config["locations"]["output-dir"], "work/refGenome"),
        )
    except FileExistsError:
        pass

    # Create file links
    for sample in config["SAMPLES"]:
        flist = config["SAMPLES"][sample]["files"]
        single_end = len(flist) == 1

        # TODO: trimming
        if not config["general"]["trimming"]:
            for idx, f in enumerate(flist):
                if not f.endswith(".gz"):
                    # FIXME: Future versions should handle unzipped .fq or .bz2.
                    bail("Input files must be gzipped: %s." % f)

                tag = "_trimmed" if single_end else f"_{str(idx+1)}_val_{str(idx+1)}"
                linkname = config["SAMPLES"][sample]["SampleID"] + tag + ".fq.gz"
                makelink(
                    path.join(config["locations"]["input-dir"], f),
                    path.join(
                        config["locations"]["output-dir"], "work/input/", linkname
                    ),
                )
        else:
            for idx, f in enumerate(flist):
                if not f.endswith(".gz"):
                    # FIXME: Future versions should handle unzipped .fq or .bz2.
                    bail("Input files must be gzipped: %s." % f)

                tag = "" if single_end else "_" + str(idx + 1)
                linkname = config["SAMPLES"][sample]["SampleID"] + tag + ".fq.gz"
                makelink(
                    path.join(config["locations"]["input-dir"], f),
                    path.join(
                        config["locations"]["output-dir"], "work/input/", linkname
                    ),
                )


def cleanup():
    out_dir = config["locations"]["output-dir"]
    need_to_rm_dirs = ["01_raw_QC", "02_trimming", "03_posttrimming_QC", "04_mapping",
                       "05_sorting_deduplication", "06_methyl_calls", "07_bigwig_files", "08_segmentation", "Reports", "09_differential_methylation","work"]
    need_to_rm_files = ["colData.tsv"]

    for file in need_to_rm_files:
        file_path = os.path.join(out_dir, file)
        if os.path.exists(file_path):
            print(f"Looks like a previous output file exists, remove files: {file_path}")
            os.remove(file_path)

    for dir in need_to_rm_dirs:
        path = os.path.join(out_dir, dir)
        if os.path.exists(path):
            print(f"Looks like a previous output file exists, remove dir: {path}")
            shutil.rmtree(path)
    print("Clean up done!!")


# Init?
if args.init:
    init_settings = False
    init_sample_sheet = False
    base = path.join(WORKDIR, "etc")

    if args.init == "both":
        init_settings = True
        init_sample_sheet = True

    if args.init == "settings" or init_settings:
        name = "settings.yaml"
        if path.exists("config/settings.yaml"):
            print("Refusing to overwrite existing {}.".format(name))
        else:
            with open(name, "w") as outfile:
                with open(path.join(base, name), "r") as infile:
                    for line in infile:
                        outfile.write("# " + line)
            os.system(f"mv {name} ./config")
            print("Generated {} template.".format(name))

    if args.init == "sample-sheet" or init_sample_sheet:
        name = "sample_sheet.csv"
        if path.exists("config/sample_sheet.csv"):
            print("Refusing to overwrite existing sample_sheet.csv.")
        else:
            shutil.copy(path.join(base, name + ".example"), name)
            os.system(f"mv {name} ./config")
            print("Generated {} template.".format(name))
    exit(0)

# Run snakemake!
generate_config(args.configfile, "config/sample_sheet.csv", "config/settings.yaml")

config = json.load(open(args.configfile, "r"))

command = [
    "snakemake",
    "-p",
    "--use-conda",
    # "--snakefile={}/snakefile.py".format(config['locations']['pkglibexecdir']),
    "--configfile={}".format(args.configfile),
    "--directory={}".format(config["locations"]["output-dir"]),
    "--jobs={}".format(config["execution"]["jobs"]),
    "--rerun-incomplete",
    "--resources",
    "mem_mb={}".format(int(config["execution"]["memory"])),
]


if args.clean:
    cleanup()
    exit(0)


if args.graph:
    command.append("--dag")
    with open(args.graph, "w") as outfile:
        p1 = subprocess.Popen(command, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["dot", "-Tpdf"], stdin=p1.stdout, stdout=outfile)
        p2.communicate()
else:
    prepare_links()
    if args.conda_create_envs_only:
        command.append("--conda-create-envs-only")
    if args.force:
        command.append("--forceall")
    if args.dry_run:
        command.append("--dryrun")
    if args.reason:
        command.append("--reason")
    if args.unlock:
        command.append("--unlock")
    if args.verbose:
        command.append("--verbose")
    if args.target and "help" in args.target:
        command.append("help")
    subprocess.run(command)
