#Make a fasta file that includes transcriptional strand
#Writes the fasta with +|-|*|N for transribed on positive, negative, both or neither strand. 
from __future__ import division, print_function
import argparse, re, sys, gzip, textwrap
from pyfaidx import Fasta
from collections import defaultdict
import numpy as np 
import pdb

MAP=np.array(["N", "+", "-", "*"])

################################################################################

def parse_options():
    """
    argparse
    """
    parser=argparse.ArgumentParser()
    parser.add_argument('-r', '--ref', type=str, default="", help="reference fasta")
    parser.add_argument('-s', '--strand', type=str, default="", help="strand file - ucsc genes track")
    return parser.parse_args()

################################################################################

def read_strand_file(strand_file):
    """
    Read the strand file information. This should be the UCSC knowngenes track info
    """

    if strand_file[-3:]==".gz":
        sfile=gzip.open(strand_file, "r")
    else:
        sfile=open(strand_file, "r")

    header=sfile.next()
    bits=header.split()
    if not (bits[1]=="chrom" and bits[2]=="strand" and bits[3]=="txStart" and bits[4]=="txEnd"):
        raise Exception("colums 2-5 should be chrom,strand,txStart,txEnd")

    transcripts=defaultdict(list)
    for row in sfile:
        bits=row.split()
        transcripts[bits[1]].append(bits[2:5])
    
    return transcripts
    
###############################################################################

def main(options):
    transcripts=read_strand_file(options.strand)
    ref=Fasta(options.ref)
    for chrom in ref.keys():
        print(">"+chrom)
        plus=np.array([False]*len(ref[chrom]))
        minus=np.array([False]*len(ref[chrom]))

        for transcript in transcripts[chrom]:
            if transccipt[0]=="+":
                plus[transcript[1]:transcript[2]]=True
            elif transcript[0]=="-":
                minus[transcript[1]:transcript[2]]=True

        chrom_tx_strand = "".join(MAP[1*plus+2*minus])
        print(textwrap.fill(chrom_tx_strand,40))

################################################################################

options=parse_options()
main(options)
