# Make a fasta file with deciles of recombination rate

from __future__ import division, print_function
import argparse, re, sys, gzip
from pyfaidx import Fasta
from collections import defaultdict
from itertools import chain
from scipy import interpolate
import numpy as np 
import pdb

BASES=np.array(["a", "c", "g", "t", "A", "C", "G", "T"])
AUTOSOMES=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
           "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"]


################################################################################

def parse_options():
    """
    argparse
    """
    parser=argparse.ArgumentParser()
    parser.add_argument('-r', '--ref', type=str, default="", help="reference fasta")
    parser.add_argument('-c', '--rec', type=str, default="", help="path to recombination maps /genetic_map_GRCh37_chr*.txt.gz")
    parser.add_argument('-b', '--block', type=int, default=10000, help="block size")
    return parser.parse_args()

################################################################################

class recombinator(object):
    """
    Class to interpolate the genetic distance between two positions, given am IMPUTE stype recombination map.
    This is c_j d_j from Li and Stephens. 
    """
    
    def __init__(self, recombination_map):
        
        recombination_file=None
        if recombination_map[-3:]==".gz":
            recombination_file=gzip.open(recombination_map, "r")
        else:
            recombination_file=open(recombination_map, "r")
        
        self.position=[]
        self.rate=[]
        self.dist=[]

        header=recombination_file.next()
        header_map=self.parse_header(header)
        
        for line in recombination_file:
            data=line[:-1].split()
            self.position.append(int(data[header_map["POS"]]))
            self.rate.append(float(data[header_map["RATE"]]))
            self.dist.append(float(data[header_map["MAP"]]))
        
        recombination_file.close()

        self.max_pos = max(self.position)
        self.min_pos = min(self.position)

        self.fitter = interpolate.UnivariateSpline(self.position, self.dist, k=1, s=0) # Linear interpolation

    def parse_header(self, header_line):
        """
        Try and parse the header line and figure out what are the position, rate and dist columns
        """
        bits=header_line.split()
        
        header_map={"POS":None, "RATE":None, "MAP":None}
        
        for what in header_map.keys():
            where=[y.upper()[:len(what)]==what for y in bits]
            if sum(where)!=1:
                raise Exception("Check that map file has headers \"Pos\", \"Rate\" and \"Map\" (in CM)")
            else:
                header_map[what]=[i for i in range(len(where)) if where[i]][0]

        return header_map
                
    def distance(self, position_1, position_2):
        """
        Return the genetic distance in cm between two points
        """
        map_pos = self.fitter([position_1, position_2])
        return map_pos[1] - map_pos[0]

#########################################################################################################

def main(options):
    
    ref=Fasta(options.ref)

    blocks={}
    for chrom in AUTOSOMES:
        blocks[chrom]=[]
        chrlen=len(ref[chrom])
        rec_map=recombinator(options.rec+"/genetic_map_GRCh37_chr"+chrom+".txt.gz")
        i=0
        while options.block *(i+1) < chrlen:
            r_rate=rec_map.distance(i, options.block*(i+1))
            blocks[chrom].append(r_rate)
        
        #Last chunk
        last_chunk_start=options.block *(i+1)
        last_chunk_end=chrlen
        r_rate=rec_map.distance(last_chunk_start, last_chunk_end)
        r_rate=r_rate*options.block/(last_chunk_start-last_chunk_end)
        blocks[chrom].append(r_rate)
        blocks[chrom]=np.array(blocks[chrom])
        
    all_rates=np.array(chain(*blocks.values()))
    quantiles=np.percentile(all_rates, 10*np.array(range(1,11)))
    quantile_blocks={}
    for k,v in blocks:
        quantile_blocks[k]=np.searchsorted(quantiles, v)
        
    #Write out result
    for chrom in AUTOSOMES:
        print(">"+chrom)
        seq=ref[chrom].seq
        quantile_chr=[str(x)*options.block for x in quantile_blocks[chrom]]
        quantile_chr=quantile_chr[:len(seq)]
        out_qsec="".join([q if x in BASES else "N" for x,q in zip(seq, quantile_chr)])
        print(out_qsec)
        
        
    
################################################################################

if __name__=="main":
    options=parse_options()
    main(options)
