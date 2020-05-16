#!/usr/bin/env python3

__title__        = 'TR_ref_gen'
__version__      = '0.2.0'
__description__  = "Generate a VaNTA Tandem Repeat reference library from a set of reference genes"
__author__       = 'Arthur V. Morris'
__institution__  = "IBERS - Aberystwyth University"
__author_email__ = "arthurvmorris@gmail.com | arm21@aber.ac.uk"

__doc__ = " %s v%s by %s\n - %s - \nContact: %s\n%s" % (__title__,
                                                     __version__,
                                                     __author__,
                                                     __description__,
                                                     __author_email__,
                                                     __institution__)

import sys
import os
import argparse
import subprocess
import shutil
import re
import collections as c
from VaNTAheader import *
from pipes import quote
from Bio import SeqIO, Seq

from VaNTAheader import *

subprocessID = "TR-GEN"

def make_fasta(gene_dict, args, gene_fa_out):

    with open(gene_fa_out, "w") as out_handle:
        for rec in SeqIO.parse(args.fasta, "fasta"):
            for gene, start, end, strand in gene_dict[rec.id]:
                seq = str(rec.seq)[start-(args.flank+1):end+args.flank]
                if strand == "-":
                    seq = Seq.reverse_complement(seq)

                out_handle.write(f">{gene} | location={rec.id}:{start}-{end} | strand={strand}\n{seq}\n")


def run_trf(args, gene_multi, trf_out_file="./VaNTA_trf.dat"):

    trf_runline = "{trf_path} {gene_multi} {trf_args} >> {trf_out} 2>&1".format(trf_path=args.TRF, \
                                                                                gene_multi=gene_multi, \
                                                                                trf_args=args.trf, \
                                                                                trf_out=trf_out_file)

    # trf_out_file = "{output_dir}/{basename}.{args}.dat".format(output_dir=output_dir, basename=os.path.basename(gene_multi), args=".".join(trf_args[:7]))

    origWD = os.getcwd()

    os.chdir(os.path.abspath(args.output))

    if args.suppress is False:
        vprint(subprocessID,
         f"Running TRF with runline: {trf_runline}\n",
         "prYellow")
        subprocess.call(trf_runline, shell=True)
    else:
        vprint(subprocessID,
         f"Assuming trf .dat file is {os.path.abspath(trf_out_file)}",
         "prYellow")

    os.chdir(os.path.abspath(origWD))

    return trf_out_file

def read_trf_dat(output_file, args, trf_out_file="./VaNTA_trf.dat"):

    """Parses trf output, where:
    start, stop = Indices of the repeat relative to the start of the sequence.
    PS = Period size of the repeat.
    CN = Number of copies aligned with the consensus pattern.
    SC = Size of consensus pattern (may differ slightly from the period size).
    PM = Percent of matches between adjacent copies overall.
    PI = Percent of indels between adjacent copies overall.
    AS = Alignment score.
    P_{A/C/G/T} = Percent composition for each of the four nucleotides.
    EM = Entropy measure based on percent composition.
    """

    tr_gff_out = open(output_file, "w")

    with open(trf_out_file, "r") as tof:

        g_TR_recs = tof.read().split("@")

        for g_rec in g_TR_recs:
            g_rec_tmp = g_rec.split("\n")
            header = g_rec_tmp[0]
            if header != "":
                gene_id, _, _, _, chr_id, gene_interval, _, _, _, strand = re.split("[||:|=| |\n]", header)
                gene_start, gene_stop = gene_interval.split("-")
                for tr_rec in (g for g in g_rec_tmp[1:] if g != ""):
                    tr_rec_tmp = tr_rec.split(" ")

                    start, \
                    end, \
                    PS, \
                    CN, \
                    SC, \
                    PM, \
                    PI, \
                    AS, \
                    P_A, \
                    P_C, \
                    P_G, \
                    P_T, \
                    EM, \
                    period_seq, \
                    rep_seq, \
                    f1_seq, \
                    f2_seq = tr_rec_tmp

                    description = "ID={gene_id}.P.{tr_start}-{tr_stop};Parent={gene_id};description=tandem repeat;Period_seq={period_seq};Period_size={PS};Repeat_size={RS};Copy_number={CN};TR_seq={rep_seq};f1_seq={f1_seq};f2_seq={f2_seq}".format(\
                    gene_id=gene_id,
                    tr_start=start,
                    tr_stop=end,
                    period_seq=period_seq,
                    PS=PS,
                    RS=int(end)-int(start),
                    CN=CN,
                    rep_seq=rep_seq,
                    f1_seq=f1_seq[-args.flank:],
                    f2_seq=f2_seq[:args.flank])

                    if len(f1_seq) < args.flank or len(f2_seq) <args.flank:
                        ## checks that entire flanks are present, skips if not
                        continue
                    else:
                        gff_line = "{chr}\tVaNTA\tTR\t{tr_start}\t{tr_stop}\t.\t{strand}\t0\t{description}\n".format(\
                        chr=chr_id, tr_start=int(gene_start)+int(start), tr_stop=int(gene_start)+int(end), strand=strand, description=description)

                        tr_gff_out.write(gff_line)
                        # db_feature = gffutils.feature.feature_from_line(gff_line)
                        # gffutils.FeatureDB.update(db_feature)

    tr_gff_out.close()

def gen_argparser(argv):
    parser = argparse.ArgumentParser(description=__doc__, usage=f"{__title__}.py fasta trf_path <options>")

    parser.add_argument('script_path', action='store', help=argparse.SUPPRESS)
    parser.add_argument('fasta', action='store',
                        help='Genome sequence in fasta format.')
    parser.add_argument('gff', action='store',
                        help='File containing annotations in gff format.')
    parser.add_argument('TRF', action='store',
                        help='Path to tandem repeat finder software.')

    parser.add_argument('-a', type=str, action='store', default="all",
                        help='Specify the feature type to extract (exact string match) [Optional]. Default = extract all.')
    parser.add_argument('-o', '--output', type=str, default='./', action='store',
                        help='Path to output directory.')
    parser.add_argument('-f', '--flank', type=int, default=40, action='store',
                        help='Size of flanking regions, up to 50 nucleotides. Default=40.')
    parser.add_argument('-p', '--prefix', type=str, default='TR_ref', action='store',
                        help='gff output prefix. <PREFIX>.gff.')
    parser.add_argument('-s', '--suppress', action='store_true', default=False,
                        help='Suppress TRF run. This will look for the .dat trf output file in the provided output directory.')
    parser.add_argument('-t', '--trf', nargs=7, default='2 5 5 80 10 50 15',
                        help='TRF argument line. See the TRF manual for more details. Default = 2 5 5 80 10 50 15. DO NOT include optional arguments, -h and -ngs included.')

    args = parser.parse_args(argv)

    if args.flank > 50 or args.flank < 1:
        vprint(subprocessID,
        "Flank length specified out of range. Please specify a flank size between 1 and 50 nucleotides.",
        "prRED",
        sys.stderr)
        exit(4)

    args.trf += ' -h -ngs'   ## add compulsory trf arguments to trf argline

    return args

def main(argv):

    args = gen_argparser(argv)

    gene_fa_out = os.path.join(args.output, f"{args.prefix}_{args.a}.fa")
    tr_gff_out = os.path.join(args.output, f"{args.prefix}_TR.gff")    ## Path to TR feature file in gff format

    ### EXTRACT RECORDS FROM GFF AND GENERATE SEQRECORD ###
    gene_dict = parse_gff(args)
    make_fasta(gene_dict, args, gene_fa_out)

    ### RUN TANDEM REPEATS FINDER USING SEQRECORD FEATURE LIBRARY AND TRF ARGLINE ###
    trf_out_file = run_trf(args, gene_fa_out)

    ### READ TRF OUTPUT FILE AND GENERATE A TANDEM REPEAT ANNOTATION FILE IN GFF FORMAT ###
    read_trf_dat(tr_gff_out, args)

genTRrefMAIN = main

if __name__ == "__main__":
    main(sys.argv)
