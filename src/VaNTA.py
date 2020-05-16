#!/usr/bin/env python

"""
888     888         888b    888 88888888888     d8888
888     888         8888b   888     888        d88888
888     888         88888b  888     888       d88P888
Y88b   d88P 8888b.  888Y88b 888     888      d88P 888
 Y88b d88P     "88b 888 Y88b888     888     d88P  888
  Y88o88P  .d888888 888  Y88888     888    d88P   888
   Y888P   888  888 888   Y8888     888   d8888888888
    Y8P    "Y888888 888    Y888     888  d88P     888

VNTR discovery pipeline v1.0
Arthur V. Morris 2020
arthurvmorris@gmail.com | arm21@aber.ac.uk
"""

import sys
import argparse, argcomplete
import subprocess
import re
import random, string
from time import localtime, strftime
from os import mkdir, path

from VaNTAheader import *
from genTRref import *
from VaNTAalign import *

subprocessID = "VaNTA"

class Vargs:
    def __init__(self, config_file):
        argv = read_config_file(config_file)

        self.script_dir = argv[0]
        self.wd = is_dir(argv[1])
        self.ref_gff = is_file(argv[2])
        self.ref_fa = is_file(argv[3])
        self.query_dir = is_dir(argv[4])
        self.feature_type = argv[5]
        self.trf_path = argv[6]
        self.flank_size = int(argv[7])
        self.VNTR_out_prefix = argv[8]
        self.water_path = argv[9]

def read_config_file(config_file):
    argv = []
    config_file = open(config_file, "r")

    if config_file.readline() != "##VaNTA_config_file##\n":
        vprint(subprocessID,
         "Control file header not recognised! Please check it is correctly formatted. Exiting...",
          "prRed",
          sys.stderr)
        sys.exit(2)

    lines = config_file.readlines()

    for line in filter(lambda w: w.startswith("~#"), lines):
        line = re.split("[ \n]", line)[:-1]
        line = line[1:]

        argv.append(line[0])

    return argv

def gen_file_structure(vargs, args):
    """ Generates the file structure to deposit VaNTA run files.
    This file structure is as follows:

         working_directory
                 |
            VaNTA.runkey
                 |
        -------------------
        |        |        |
     trf_out  aln_out multifastas

    Generates files with a randomly generated run key, which is (hopefully) unique to this run.

    """
    ## generates a unique run code to name directories, a bit hacky...
    key = genKey()
    vprint(subprocessID,
     f"Run key for this run: {key}\n",
      "prYellow")

    vargs.runWD = path.join(vargs.wd, f"VaNTA.{key}")
    mkdir(vargs.runWD)

    vprint(subprocessID,
     f"Generating directory structure for this VaNTA run in working directory {vargs.runWD}\n",
      "prYellow")

    vargs.trf_out = path.join(vargs.runWD, "trf_out")
    if not path.isdir(vargs.trf_out):
        mkdir(vargs.trf_out)   ## make trf out dir
    else:
        vprint(subprocessID,
         f"{vargs.trf_out} exists! Aborting to prevent accidental data loss...\n",
          "prRed",
          sys.stderr)
        exit(3)

    vargs.aln_out = path.join(vargs.runWD, "aln_out")
    if not path.isdir(vargs.aln_out):
        mkdir(vargs.aln_out)   ## make directory to contain alignment files
    else:
        vprint(subprocessID,
         f"{vargs.aln_out} exists! Aborting to prevent accidental data loss...\n",
          "prRed",
          sys.stderr)
        exit(3)

    vargs.multi_out = path.join(vargs.runWD, "multifastas")
    if not path.isdir(vargs.multi_out):
        mkdir(vargs.multi_out) ## make directory to contain multifastas to parse to water aligner
    else:
        vprint(subprocessID,
         f"{vargs.multi_out} exists! Aborting to prevent accidental data loss...\n",
          "prRed",
          sys.stderr)
        exit(3)

def gen_config_file():
    """ Generates a blank control file to parse to VaNTA.py main script.
    """
    out_handle = open("./VaNTA_args.txt", "w")

    args_string = """##VaNTA_config_file##

~#script_dir: /home/user/VaNTA/script/
~#working_directory: /home/user/VaNTA
~#ref_gff: /home/user/VaNTA/ref/ref.gff
~#ref_fa: /home/user/VaNTA/ref/ref.fa
~#query_dir: /home/user/VaNTA/query_lib/
~#feature_type: gene
~#TRF_bin: /home/user/VaNTA/script/trf409.linux64
~#flank_size: 40
~#VNTR_multifasta_out_prefix: ref.VNTR
~#EMBOSS_water_path: /home/user/VaNTA/script/water

##ENDFILE##
    """

    print("Generating blank control file: ./VaNTA_args.txt")
    print(args_string, file=out_handle)
    sys.exit(0)

def genKey(): return ''.join(random.choices(string.ascii_letters + string.digits, k=6))

def runPipeline(vargs, args):
    """ Runs the full VaNTA pipeline using arguments provided in the config file
    """

    genTRrefMAIN([
        path.join(vargs.script_dir, "genTRref.py"),
        vargs.ref_fa,
        vargs.ref_gff,
        vargs.trf_path,
        "-a", vargs.feature_type,
        "-o", vargs.runWD,
        "-f", str(vargs.flank_size),
        "-p", vargs.VNTR_out_prefix
        ])

    VaNTAalignMAIN([
        path.join(vargs.script_dir, "VaNTAalign.py"),
        path.join(vargs.runWD, vargs.VNTR_out_prefix + "_TR.gff"),
        vargs.ref_fa,
        vargs.query_dir,
        "-o", vargs.aln_out,
        "-m", vargs.multi_out,
        "-c", vargs.runWD
        ])

def gen_argparser(argv):
    """ Parses arguments for VaNTA control script
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('script_path', action='store', help=argparse.SUPPRESS)
    parser.add_argument('-c', '--control_file', action="store", help='A control file containing files and parameters.')
    parser.add_argument('-g', '--gen_control', action='store_true', help='Generate a blank control file.')

    argcomplete.autocomplete(parser)
    args = parser.parse_args(argv)

    if args.gen_control and args.control_file != None:
        gen_config_file()
    elif args.control_file != None:
        vargs = Vargs(args.control_file)
        return vargs, args
    else:
        print("Please provide a VaNTA config file with the -c flag, or generate a blank one with -g.", file=sys.stderr)
        sys.exit(1)

def main(argv):

    vargs, args = gen_argparser(argv)     ## parse arguments and config file and generate an instance of Vargs containing args found in the config file
    gen_file_structure(vargs, args)       ## generates the file structure for this run of VaNTA
    runPipeline(vargs, args)


if __name__ == "__main__":
    print(__doc__)
    main(sys.argv)
