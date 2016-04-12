# Read a fasta file, and remove a motif
# Side effect of making everything uppercase
from __future__ import division, print_function
import argparse, re, sys
from pyfaidx import Fasta

motif="CC[A-Z]CC[A-Z]T[A-Z][A-Z]CC[A-Z]C"
motif_length=13
################################################################################

def parse_options():
    """
    argparse
    """
    parser=argparse.ArgumentParser()
    parser.add_argument('-r', '--ref', type=str, default="", help="reference fasta")
    return parser.parse_args()

################################################################################

def main(options):
    """
    Iterate and remove motif (by setting to N)
    """
    ref=Fasta(options.ref)
    reg=re.compile(motif)
    for chrom in ref.keys():
        print(">"+chrom)
        new_seq=reg.sub("N"*motif_length, ref[chrom][:].seq.upper())
        print(new_seq)
        
################################################################################

options=parse_options()
main(options)
