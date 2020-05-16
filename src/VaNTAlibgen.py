#!/usr/bin/env python3

import sys
import argparse
import re
from VaNTAheader import *
from os import path
from Bio import SeqIO, Seq


__title__        = 'VaNTA_library_generator'
__version__      = '0.1.2'
__description__  = "Generate VaNTA CDS libraries from .vff files"
__author__       = 'Arthur V. Morris'
__institution__  = "IBERS - Aberystwyth University"
__author_email__ = "arthurvmorris@gmail.com | arm21@aber.ac.uk"

__doc__ = " %s v%s by %s\n - %s - \nContact: %s\n%s" % (__title__,
                                                     __version__,
                                                     __author__,
                                                     __description__,
                                                     __author_email__,
                                                     __institution__)

def make_fasta(gene_dict, args):

    with open(args.out_file, "w") as out_handle:
        for rec in SeqIO.parse(args.fasta, "fasta"):
            for gene, start, end, strand in gene_dict[rec.id]:
                seq = str(rec.seq)[start-args.f-1:end+args.f]
                if strand == "-":
                    seq = Seq.reverse_complement(seq)

                out_handle.write(f">{gene} | location={rec.id}:{start-args.f-1}-{end+args.f} including {args.f}n flank | strand={strand}\n{seq}\n")

def execute(argv):

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('script_path', action='store', help=argparse.SUPPRESS)
    parser.add_argument('gff', action='store', help='The annotation file in gff format')
    parser.add_argument('fasta', action='store', help='The sequence file in FASTA format.')
    parser.add_argument('out_file', action='store', help='Output file.')

    parser.add_argument('-a', type=str, action='store', default="all", help='Specify the feature type to extract (exact string match) [Optional]. Default = extract all.')
    parser.add_argument('-f', type=int, action='store', default=40, help='Size of the flanking regions to include when extracting gene sequences. Default = 40.')

    args = parser.parse_args(argv)

    sys.stdout.write(f"Extracting {args.a} features to {args.out_file} with flank size {args.f}\n")

    gene_dict = parse_gff(args)
    make_fasta(gene_dict, args)

execute(sys.argv)
# execute(["foo", "/home/amorris/BioInf/VaNTA2_WD/CryptoDB-46_CparvumIowaII.gff", "/home/amorris/BioInf/VaNTA2_WD/CryptoDB-46_CparvumIowaII_Genome.fasta", "/home/amorris/BioInf/VaNTA2_WD/IowaII_gene.fa", "-a", "gene"])
