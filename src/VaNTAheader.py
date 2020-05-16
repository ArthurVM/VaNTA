"""
A header file for VaNTA scripts containing generic functions
"""

import sys
import argparse
import re
import collections as c
from time import localtime, strftime
from Bio import SeqIO
from os import path

colfuncs = {}
colfunc = lambda f: colfuncs.setdefault(f.__name__, f)

@colfunc
def prRed(sp):     return f"\033[91m{sp}\033[00m"
@colfunc
def prGreen(sp):   return f"\033[92m{sp}\033[00m"
@colfunc
def prYellow(sp):  return f"\033[93m{sp}\033[00m"

def is_file(filename):
    """ Checks if a path is a file """

    if not path.isfile(filename):
        time=strftime("%H:%M:%S", localtime())
        print("{time} {subprocess} :: No file found at {f}".format(
                time=time,
                subprocess=prRed("ERROR"),
                f=filename), end="\n", file=sys.stderr, flush=True)
        exit(3)
    else:
        return path.abspath(path.realpath(path.expanduser(filename)))

def is_dir(dirname):
    """ Checks if a path is a directory """

    if not path.isdir(dirname):
        time=strftime("%H:%M:%S", localtime())
        print("{time} {subprocess} :: No directory found at {d}".format(
                time=time,
                subprocess=prRed("ERROR"),
                d=dirname), end="\n", file=sys.stderr, flush=True)
        exit(3)
    else:
        return path.abspath(path.realpath(path.expanduser(dirname)))

def file_id(file_path):
    return path.splitext(path.basename(file_path))[0]

def vprint(subprocess, info_text, colour, f=sys.stdout, end="\n"):
    """ controls process output
    """

    time=strftime("%H:%M:%S", localtime())

    print("{time} {subprocess} :: {info_text}".format(
            time=time,
            subprocess=colfuncs[colour](subprocess),
            info_text=info_text),
            end=end, file=f, flush=True)

def parse_gff(args):

    gene_dict = c.defaultdict(list)

    with open(args.gff, "r") as in_handle:
        for line in in_handle:
            if line.startswith("##"):
                continue

            sline = line.split("\t")

            if sline[2] == args.a or args.a == "all":
                chr = sline[0]

                desc = sline[-1].split(";")
                for d in desc:
                    if d.startswith("ID"):
                        id = re.split("[=|\n]", d)[1]

                gene_dict[chr].append([id, int(sline[3]), int(sline[4]), sline[6]])

    return gene_dict
