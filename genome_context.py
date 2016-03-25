#Get the distribution of trinucleotide contexts from a reference genome.

from __future__ import division, print_function
import argparse, sys, itertools, re
from pyfaidx import Fasta

BASES=["A","C","G","T"]
BASE_MAP={"A":("A","A"), "C":("C","C"), "G":("G", "G"), "T":("T","T"), "M":("A", "C"), "R":("A","G"), "W":("A","T"), "S":("C","G"), "Y":("C", "T"), "K":("G","T")}
CHRS=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"]
################################################################################

def parse_options():
    """
    argparse
    """
    parser=argparse.ArgumentParser()
    parser.add_argument('-r', '--ref', type=str, default="", help="reference fasta")
    parser.add_argument('-k', '--kmer', dest="k", type=int, default=3, help="kmer to consider")
    parser.add_argument('-d', '--hetfa', action='store_true', help="Input diploid fasta (hetfa)")
    parser.set_defaults(hetfa=False)

    return parser.parse_args()

##########################################################################################################

def make_dict(k):
    '''
    Make the structure to store the results - a dict with keys given by ('XY...Z', 0) which
    represents the count of the "XY...Z" tuple. 
    '''    
    keys=["".join(x) for x in itertools.product(*[BASES]*k)]
    d={}.fromkeys(keys, 0)                         
    return d


################################################################################

def main(options):
    '''
    Scan through the reference genome and output the counts of trinucleotide contexts
    '''
    ref=Fasta(options.ref)
    d=make_dict(options.k)
    allowed_chars=reg=re.compile( "^[ACGTMRWSYK]+$" )
    for chrom in CHRS:
        print("chr"+chrom, file=sys.stderr)
        chr_length=len(ref[chrom])
        for i in xrange(chr_length-options.k):
            tnc=ref[chrom][i:(i+options.k)].seq.upper()
            if not options.hetfa and tnc in d:
                d[tnc]+=1
            elif options.hetfa and allowed_chars.match(seq):
                c1="".join((BASE_MAP[x][0] for x in tnc))
                c2="".join((BASE_MAP[x][1] for x in tnc))
                d[tnc]+=1
                
    for k,v in d.iteritems():
        print(k+"\t"+str(v))
    
################################################################################

options=parse_options()
main(options)
