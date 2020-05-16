"""
Tandem Repeat class.

The tandem repeat (TR) class contains
"""

import sys
import os
import argparse
import subprocess
import datetime
import re
import numpy as np
from VaNTAheader import *
from Bio import SeqIO, Seq
import collections as c

class TR(object):
    """docstring for tandem_repeat class."""

    def __init__(self, gff_line, fasta_dict, flank):
        gl = gff_line
        desc = parse_desc(gl[-1])

        self.seq = desc["seq"]
        self.flank1 = desc["f1"]
        self.flank2 = desc["f2"]
        self.chr = gl[0]
        self.start = int(gl[3])
        self.end = int(gl[4])
        self.strand = gl[6]
        self.id = desc["id"]
        self.parent = desc["parent"]
        self.period_seq = desc["period_seq"]
        self.size = int(desc["size"])
        self.copy_number = float(desc["cn"])

    def get_homologs(self, seqrec_array, out_dir=False):
        """ Generates a dictionary containing homologous sequences to align from the seqrec_array using self.parent.
        Each element within seqrec_array should be a biopython sequence record generated from a query fasta file.
        """
        homdict = c.defaultdict(str)
        missing = 0
        for id, recs in seqrec_array.items():
            if self.parent in recs:
                homdict[id]=recs[self.parent]
            else:
                missing+=1
                continue
        return homdict, missing
        # return { id : recs[self.parent] if self.parent in recs else None for id, recs in seqrec_array.items() }

    def get_homologs_fasta(self, seqrec_array, out_dir):
        """ Generates 2 fasta files:
                TR reference seqeuence
                Query homolog multifasta
            and outputs to out_dir.
        """
        missing = 0
        rfa_out = path.join(out_dir, f"{self.id}.fa")
        rfa_F1_out = path.join(out_dir, f"{self.id}.F1.fa")
        rfa_F2_out = path.join(out_dir, f"{self.id}.F2.fa")
        qfa_out = path.join(out_dir, f"{self.parent}.fa")
        with open(rfa_out, "w") as ref_handle, open(qfa_out, "w") as q_handle:
            ref_handle.write(f">{self.id}\n{self.flank1}{self.seq}{self.flank2}\n>{self.id}_flank1\n{self.flank1}\n>{self.id}_flank2\n{self.flank2}\n")
            for id, recs in seqrec_array.items():
                if self.parent in recs:
                    q_handle.write(f">{id}\n{str(recs[self.parent].seq)}\n")
                else:
                    missing+=1
                    continue
        return rfa_out, qfa_out, missing

    def parse_aln(self, aln_dict, qlib_n):
        """ Takes an intance of TR and a dictionary containing alignment data and updates the instance.
        aln_dict: a dictionary containing alignment data structured as
                    {   full_seq : { Q0 : match, ..., Qn : match },
                        flank1   : { Q0 : match, ..., Qn : match },
                        flank2   : { Q0 : match, ..., Qn : match }    }

        This is used to generate the attributes:
            allele  : number of alleles in dataset
            $_avg   : mean percentage match     where $=tr,f1,f2
            count   : number of queries in which the TR is found
            V_score : the V-score for this TR locus

        The scores where produced using this method:

                SNP penality = 5*(similarity score + gaps)

        This is to penalise sequences where the variability is due to SNP's rather than gaps, which
        are assumed to be due to VNTR copy number variation.

        VNTR variability score (vv_score) is then calculated using the following equation:

                B + V + f - l

                where: B = base score
                       V = TR variability score
                       f = flank penalty
                       l = TR length
        """
        fullseq_array = [match for id, match in aln_dict["full_seq"].items()]
        f1_array = [match for id, match in aln_dict["flank1"].items()]
        f2_array = [match for id, match in aln_dict["flank2"].items()]

        TR_allele_n = len(set(fullseq_array))

        V_score = 100*(TR_allele_n) + (np.mean(f1_array)+np.mean(f2_array)) - (10*len(self.period_seq)) - len(self.seq)

        # print(100*(TR_allele_n), (np.mean(f1_array)+np.mean(f2_array)), (10*len(self.period_seq)), len(self.seq))

        self.allele = TR_allele_n
        self.tr_avg = np.mean(fullseq_array)
        self.f1_avg = np.mean(f1_array)
        self.f2_avg = np.mean(f2_array)
        self.count = len(fullseq_array)
        self.V_score = V_score
        self.libn = len(fullseq_array)

    def tsv_line(self):
        """ takes an instance of TR and generates a tsv line.
        """
        if hasattr(self, "V_score"):
            return f"{self.id}\t{self.period_seq}\t{len(self.seq)}\t{self.libn}\t{self.f1_avg}\t{self.f2_avg}\t{self.allele}\t{self.V_score}"
        else:
            return None

def parse_desc(desc):
    desc_dict = {}
    tmp = re.split("[=|;]|\n", desc)
    desc_dict["id"] = tmp[tmp.index("ID")+1]
    desc_dict["parent"] = tmp[tmp.index("Parent")+1]
    desc_dict["period_seq"] = tmp[tmp.index("Period_seq")+1]
    desc_dict["size"] = tmp[tmp.index("Repeat_size")+1]
    desc_dict["cn"] = tmp[tmp.index("Copy_number")+1]
    desc_dict["seq"] = tmp[tmp.index("TR_seq")+1]
    desc_dict["f1"] = tmp[tmp.index("f1_seq")+1]
    desc_dict["f2"] = tmp[tmp.index("f2_seq")+1]

    return desc_dict

def parse_TR_gff(args):
    gff_frame = []
    with open(args.TR_gff) as fin:
        for line in fin.readlines():
            tr_line = line.split("\t")
            gff_frame.append(
                TR( tr_line,
                    SeqIO.to_dict(SeqIO.parse(args.ref_fasta, "fasta")),
                    args.f )
                )
    return gff_frame
