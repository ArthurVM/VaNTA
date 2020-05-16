#!/usr/bin/env python3

"""
VaNTAalign aligns reference TR regions to homologs within the CDS library array.
"""

# TODO: EMBOSS water reader class

import sys
import argparse
import re
import numpy as np
from subprocess import Popen, PIPE
from time import localtime, strftime
from os import path, listdir, walk
from VaNTAheader import *
from TRframe import *
from Bio.Emboss.Applications import WaterCommandline, NeedleCommandline
from Bio import SeqIO, Align, Seq

subprocessID = "ALIGN"

def tr_aligner(TR_frame, args, m=10.0, go=10.0, ge=0.5):
    """ Controls  the alignment of TR regions and flanks
    """
    vprint(subprocessID,
     f"Reading gene libraries from {args.CDS_libs}",
     "prYellow",
     f=sys.stdout,
     end="...")

    seqrec_array = { file_id(CDS_lib):SeqIO.index(path.join(args.CDS_libs, CDS_lib), "fasta")\
                for CDS_lib in listdir(args.CDS_libs) }

    print("Done.\n", flush=True)

    if args.w is False:
        missing = bio_aligner(TR_frame, seqrec_array, m, go, ge, args)

    else:
        missing = water_aligner(TR_frame, seqrec_array, m, go, ge, args)

def water_aligner(TR_frame, seqrec_array, m, go, ge, args):
    """ Performs TR alignment using the provided EMBOSS-water aligner executable.
    TR_frame: A data frame containing TR instances
    seqrec_array: An array containing indexed seqrecord instances from the query feature array
    m: match score
    go: gap open penalty
    ge: gap extension penalty
    min_match: The minimum percentage similarity to accept the alignment, otherwise realign with reverse complement or remove
    """
    tr_count = len(TR_frame)

    missing_features = 0    ## counter for instances of missing features in the query lib

    water_log = open("./TR_aln.log", "w")   ## file to dump water subprocess output

    vprint(subprocessID,
     "Starting alignments...",
     "prYellow")
    print(f"\n\t\t\tEMBOSS-water Smith-Waterman Aligner.\n\t\t\tmatch={m}\n\t\t\tgap_open={go}\n\t\t\tgap_extend={ge}\n", flush=True)

    for i, tr in enumerate(TR_frame):
        time=strftime("%H:%M:%S", localtime())
        print("\r{time} {subprocess} :: Aligning TR {i}/{tr_count}".format(
                time=time,
                subprocess=prYellow(subprocessID),
                i=i+1,
                tr_count=tr_count), end="... ", file=sys.stdout, flush=True)

        ## generate a homolog dict from the TR and track the number of missing features from the query library
        rfa_out, qfa_out, missing = tr.get_homologs_fasta(seqrec_array, args.m)
        missing_features+=missing

        aln_out = path.join(args.o, f"{tr.id}.water")
        flank1_out = path.join(args.o, f"{tr.id}_F1.water")
        flank2_out = path.join(args.o, f"{tr.id}_F2.water")

        def run_alignment(water_cline, a_prefix):
            """ Runs water alignment using a Biopython water commandline object and a prefix to identify which sequence is being aligned
            """
            p = Popen(str(water_cline), stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=True)
            output, err = p.communicate()
            rc = p.returncode

            if rc == 0:
                print(f"\nAlignment of {a_prefix}:{tr.id} exited with 0", file=water_log)
            else:
                print(f"\nAlignment of {a_prefix}:{tr.id} exited with {rc} and warning:\n{str(err)}", file=water_log)

        ## align full TR region
        water_cline = WaterCommandline(
            args.w,
            asequence=rfa_out,
            bsequence=qfa_out,
            gapopen=go,
            gapextend=ge,
            outfile=aln_out)

        run_alignment(water_cline, "FULL")

        ## align flank 1
        water_cline = WaterCommandline(
            args.w,
            asequence=f"asis:{tr.flank1}",
            bsequence=qfa_out,
            gapopen=go,
            gapextend=ge,
            outfile=flank1_out)

        run_alignment(water_cline, "flank1")

        ## align flank 2
        water_cline = WaterCommandline(
            args.w,
            asequence=f"asis:{tr.flank2}",
            bsequence=qfa_out,
            gapopen=go,
            gapextend=ge,
            outfile=flank2_out)

        run_alignment(water_cline, "flank2")

        if i > 10:
            break

    print(f"Done with {missing_features} missing seqeuences.\n", flush=True)
    water_log.close()

def bio_aligner(TR_frame, seqrec_array, m, go, ge, args, min_match=60):
    """ Performs TR alignment using the Biopython Align module.
    TR_frame: A data frame containing TR instances
    seqrec_array: An array containing indexed seqrecord instances from the query feature array
    m: match score
    go: gap open penalty
    ge: gap extension penalty
    min_match: The minimum percentage similarity to accept the alignment, otherwise discard under the assumption that it is erroneous
    """

    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = m
    aligner.open_gap_score = -go
    aligner.extend_gap_score = -ge

    tr_count = len(TR_frame)

    nan_box = []            ## box for IDs with "nan" values
    non_homo = []           ## box for non-homologous feature IDs
    missing_features = 0    ## counter for instances of missing features in the query lib

    vprint(subprocessID,
     "Starting alignments...",
     "prYellow")
    print(f"\n\t\t\tBiopython {aligner.algorithm}.\n\t\t\tmatch={m}\n\t\t\tgap_open={go}\n\t\t\tgap_extend={ge}\n", flush=True)

    for i, tr in enumerate(TR_frame):
        aln_dict = c.defaultdict(dict)  ## a dictionary for alignment data

        time=strftime("%H:%M:%S", localtime())
        print("\r{time} {subprocess} :: Aligning TR {i}/{tr_count}".format(
                time=time,
                subprocess=prYellow(subprocessID),
                i=i+1,
                tr_count=tr_count), end="... ", file=sys.stdout, flush=True)

        with open(path.join(args.o, f"{tr.id}.aln"), "w") as aln_out,\
            open(path.join(args.o, f"{tr.id}.f1.aln"), "w") as f1_aln_out,\
            open(path.join(args.o, f"{tr.id}.f2.aln"), "w") as f2_aln_out:

            print(f"TR_id={tr.id}", file=aln_out)
            full_trseq = tr.flank1+tr.seq+tr.flank2

            ## generate a homolog dict from the TR and track the number of missing features from the query library
            hom_dict, missing = tr.get_homologs(seqrec_array)
            missing_features+=missing

            for id, hom in hom_dict.items():
                if len(hom.seq) == 0:
                    missing+=1
                    continue

                aln_info, tr_perc = run_aln(aligner, full_trseq, str(hom.seq), tr, id)  ## perform TR alignment
                print(aln_info, file=aln_out)
                if tr_perc >= min_match:
                    aln_dict["full_seq"][id]=tr_perc
                else:
                    non_homo.append(tr.id)

                aln_info, f1_perc = run_aln(aligner, tr.flank1, str(hom.seq), tr, id)  ## perform flank1 alignment
                print(aln_info, file=f1_aln_out)
                if tr_perc >= min_match:
                    aln_dict["flank1"][id]=f1_perc

                aln_info, f2_perc = run_aln(aligner, tr.flank2, str(hom.seq), tr, id)  ## perform flank2 alignment
                print(aln_info, file=f2_aln_out)
                if tr_perc >= min_match:
                    aln_dict["flank2"][id]=f2_perc

        tr.parse_aln(aln_dict, len(args.CDS_libs))

        # if i == 10:
        #     break

    print(f"Done with {missing_features} missing seqeuences.\n")

    tsv_out = path.join(args.c, "VaNTA.tsv")

    with open(tsv_out, "w") as tsv_handle:

        vprint(subprocessID,
         f"Writing output to {tsv_out}",
         "prYellow",
         sys.stdout,
         "...")

        print("## VaNTA VNTR variability ##", file=tsv_handle)
        print("locus\tTR subseq\tref VNTR length\tlib_size\tmean f1 sim (%)\tmean f2 sim. (%)\tNumber of alleles\tV_score", file=tsv_handle)

        tsvline_box = []
        for tr in TR_frame:
            if tr.tsv_line() != None and "nan" not in tr.tsv_line(): tsvline_box.append(tr.tsv_line())
            else: nan_box.append(tr.id)

        stsvline_box = sorted(tsvline_box,
            key=lambda x: float(x.split("\t")[-1]), reverse=True)
        print("\n".join(stsvline_box), file=tsv_handle)

    print(f"Done.\n", flush=True)

    info_text = f"Report:\n\t\t\tTRs aligned : {len(tsvline_box)}\n\t\t\tMissing homologs : {len(non_homo)}\n\t\t\tnan rejections : {len(nan_box)}\n"
    vprint(subprocessID,
     info_text,
     "prYellow")

def run_aln(aligner, aseq, bseq, tr, id, complement=False):
    """Formats the alignment object data for output as an alignment file
    """
    alignments = aligner.align(aseq, bseq)      ## aligns the TR sequence against the query homolog
    aln = alignments[0]                         ## takes the main alignment

    r_aln, matchstr, q_aln = (i for i in re.split("[ \n]", str(aln)) if i !="")

    q_aln_interval = [np.min(aln.aligned[1]), np.max(aln.aligned[1])]
    aln_string = f"{r_aln}\n{matchstr}\n{q_aln[q_aln_interval[0]:q_aln_interval[1]]}"
    aln_info = f"{id}\nParent={tr.parent}\nLocation={tr.chr}:{tr.start}-{tr.end}({tr.strand})\nComplement adjusted={complement}"
    psl = format(aln, "psl")
    perc = 100/float(len(matchstr))*matchstr.count("|")
    aln_info = f"{aln_info}\nScore={aln.score}\nMatch={perc}\n\n{psl}\n\n{aln_string}\n"

    return aln_info, perc

def gen_argparser(argv):

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('script_path', action='store', help=argparse.SUPPRESS)

    parser.add_argument('TR_gff', type=is_file, action='store',
     help='Annotation file containing reference TR features in gff format.')
    parser.add_argument('ref_fasta', type=is_file, action='store',
     help='Reference sequence file in FASTA format.')
    parser.add_argument('CDS_libs', type=is_dir, action='store',
     help='A directory containing gene libraries.')

    parser.add_argument('-w', type=is_file, action='store', default=False,
     help='Path to the EMBOSS-water executable. If none is provided, the Biopython Align module use used to perform pairwise TR alignment.')
    parser.add_argument('-m', type=is_dir, action='store', default="./",
     help='Path to output directory for single locus multifastas generated for multiple alignment.')
    parser.add_argument('-o', type=is_dir, action='store', default="./",
     help='Path to output directory for alignment files.')
    parser.add_argument('-c', type=is_dir, action='store', default="./",
     help='Path to output directory for the .tsv file.')
    parser.add_argument('-s', action='store_true', default=False,
     help='Suppress alignments and perform only data collection [Optional]. Default=False.')
    parser.add_argument('-f', type=int, action='store', default=40,
     help='Size of the flanking regions to include when extracting gene sequences. Default = 40.')

    args = parser.parse_args(argv)

    return args

def main(argv):

    args = gen_argparser(argv)

    vprint(subprocessID,
     f"Reading in TR gff file from {args.TR_gff}",
     "prYellow",
     sys.stdout,
     "...")

    TR_frame = parse_TR_gff(args)

    print("Done.\n")

    tr_aligner(TR_frame, args)

    exit()

VaNTAalignMAIN = main
# main(sys.argv)
# main([
# "/home/amorris/BioInf/VaNTA2_WD/scripts/VaNTAalign.py",
# "/home/amorris/BioInf/VaNTA2_WD/IowaII_TR.gff",
# "/home/amorris/BioInf/VaNTA2_WD/CryptoDB-46_CparvumIowaII_Genome.fasta",
# "/home/amorris/BioInf/VaNTA2_WD/Cparv_lib",
# "-w", "/usr/bin/water",
# "-o", "/home/amorris/BioInf/VaNTA2_WD/aln_out",
# "-m", "/home/amorris/BioInf/VaNTA2_WD/multi"])
# main([
# "/home/amorris/BioInf/VaNTA2_WD/scripts/VaNTAalign.py",
# "/home/amorris/BioInf/VaNTA2_WD/IowaII_TR.gff",
# "/home/amorris/BioInf/VaNTA2_WD/CryptoDB-46_CparvumIowaII_Genome.fasta",
# "/home/amorris/BioInf/VaNTA2_WD/Cparv_lib",
# "-o", "/home/amorris/BioInf/VaNTA2_WD/aln_out",
# "-m", "/home/amorris/BioInf/VaNTA2_WD/multi"])
