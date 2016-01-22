#Get the distribution of trinucleotide contexts from a reference genome.

from __future__ import division, print_function
import argparse, sys, pdb
from pyfaidx import Fasta

BASES=["A","C","G","T"]
CHRS=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"]
################################################################################

def parse_options():
    """
    argparse
    """
    parser=argparse.ArgumentParser()
    parser.add_argument('-r', '--ref', type=str, default="", help="reference fasta")

    return parser.parse_args()

##########################################################################################################

def make_dict():
    '''
    Make the array to store the results - a dict with keys given by ('XYZ', 'W") which represents
    a trinucleotide mutation XYZ->XWZ
    '''
    
    d={}
    for a in BASES:
        for b in BASES:
            for c in BASES:
                d[a+b+c]=0
                         
    return d


################################################################################

def main(options):
    '''
    Scan through the reference genome and output the counts of trinucleotide contexts
    '''
    ref=Fasta(options.ref)
    d=make_dict()
    for chrom in CHRS:
        print("chr"+chrom, file=sys.stderr)
        chr_length=len(ref[chrom])
        for i in xrange(2, chr_length):
            tnc=ref[chrom][(i-2):(i+1)].seq
            if tnc in d:
                d[tnc]+=1
    
    for k,v in d.iteritems():
        print(k+"\t"+str(v))
    
################################################################################

options=parse_options()
main(options)
