'''
Classify the DNM mutations from GoNL
columns: ID, CHROM, POS, REF, ALT
Assume that ref is the ancestral allele (since they are de novo).
"map" remaps the IDS
@author: mathii
'''

from __future__ import division, print_function
from collections import defaultdict
from pyfaidx import Fasta
import sys, getopt, gzip, pdb

BASES=["A", "C", "G", "T"]
#Including reverse complements
SIG1=set(["TCT.T", "TCC.T", "CCC.T", "ACC.T", "AGA.A", "GGA.A", "GGG.A", "GGT.A"])

##########################################################################################################

def parse_options():
    """
    vcf: vcf input
    ref: 
    """
    options ={"dnm":None, "ref":None, "out":"results", "map":None }

    try:
        opts, args = getopt.getopt(sys.argv[1:], "d:r:o:m:",
                                    ["dnm", "ref",  "out", "map" ])

    except Exception as err:
        print( str(err), file=sys.stderr)
        sys.exit()

    for o, a in opts:
        if o in ["-d","--dnm"]:             options["dnm"] = a
        elif o in ["-r","--ref"]:           options["ref"] = a
        elif o in ["-o","--out"]:           options["out"] = a
        elif o in ["-m","--map"]:           options["map"] = a
    
    print( "found options:", file=sys.stderr)
    print( options, file=sys.stderr)

    return options, args

##########################################################################################################

def open2(f, mode="r"):
    if f.endswith(".gz"):
        return gzip.open(f, mode)
    else:
        return open(f, mode)

##########################################################################################################

def print_results(results):
    out=open(options["out"], "w")
    out.write("Mutation\t"+"\t".join(results["samples"])+"\n")
    del(results["samples"])
    for k,v in results.items():
        out.write(k[0]+"."+k[1]+"\t"+"\t".join([str(x) for x in v])+"\n")


##########################################################################################################

def read_map_file(file):
    """
    read map file - first column is value, second column is key
    """
    if not file:
        return {}

    idmap={}
    mapfile=open2(file)
    for line in mapfile:
        bits=line.split()
        idmap[bits[1]]=bits[0]

    return idmap
        
##########################################################################################################

def main(options):
    """
    Run the analysis. 
    """
    ref_error=0
        
    reference=Fasta(options["ref"])
    data=open2(options["dnm"])
    idmap=read_map_file(options["map"])
    
    total_dnm=defaultdict(int)
    sig1_dnm=defaultdict(int)

    header=data.next()
    if not all([x==y for x,y in zip(header.split()[:5], ["CHILD_ID", "CHROM", "POS", "REF", "ALT"])]):
        raise Exception("Bad headers")
        
    for line in data:
            bits=line.split()
            child_id=bits[0]
            chrom=bits[1]
            pos=int(bits[2])
            anc=bits[3]
            mut=bits[4]
            
            tnc=reference[chrom][(pos-2):(pos+1)].seq
            if tnc[1]!=anc:
                ref_error+=1
                continue
            else:
                total_dnm[child_id]+=1
                mutation=tnc+"."+mut
                if mutation in SIG1:
                    sig1_dnm[child_id]+=1
                
    print("Reference errors: "+str(ref_error))
    
    out=open(options["out"], "w")
    out.write("ID_1\tID_2\tmissing\ttotal_dnm\tsig1_dnm\tsig1_dnm_prop\n")
    out.write("0\t0\t0\tP\tP\tP\n")
    for k,v in total_dnm.items():
        mapname=idmap.get(k, k)
        out.write("\t".join([mapname, mapname, "0", str(v), str(sig1_dnm.get(k,0)), str(sig1_dnm.get(k,0)/v)])+"\n")

                         
##########################################################################################################


if __name__ == '__main__':
    options,args=parse_options()
    main(options)
