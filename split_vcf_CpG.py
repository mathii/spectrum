# Split a vcf based on whether sites are CpGs or not. 
from __future__ import division, print_function
import argparse, re, gzip
from pyfaidx import Fasta
from curses.ascii import alt

################################################################################

def parse_options():
    """
    argparse
    """
    parser=argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', type=str, default="")
    parser.add_argument('-a', '--anc', type=str, default="panTro2")    
    parser.add_argument('-r', '--ref', type=str, default="")
    return parser.parse_args()


################################################################################

def main(options):
    
    reference=Fasta(options.ref)
    
    is_gz=options.vcf[-3:]==".gz"
    
    basename=re.sub("\.vcf.+$", "", options.vcf)
    out_CpG=basename+".CpG.vcf"
    out_notCpG=basename+".nonCpG.vcf"
    if is_gz:
        out_CpG=out_CpG+".gz"
        out_notCpG=out_notCpG+".gz"
        CpG=gzip.open(out_CpG, "w")
        notCpG=gzip.open(out_notCpG, "w")
        infile=gzip.open(options.vcf, "r")
    else:
        CpG=open(out_CpG, "w")
        notCpG=open(out_notCpG, "w")
        infile=open(options.vcf, "r")

    reading_data=False
    for line in infile:        
        if not reading_data and line[:2]=="##":
            CpG.write(line)
            notCpG.write(line)
        elif not reading_data and line[0]=="#":
            CpG.write(line)
            notCpG.write(line)
            bits=line.split()
            anc_sample_idx=bits[9:].index(options.anc)
            reading_data=True
        elif reading_data:
            bits=line.split()
            bits=line.split()
            chrom=bits[0]
            pos=int(bits[1])
            ref=bits[3]
            alt=bits[4]
            gts=bits[9:]
            
            if gts[anc_sample_idx]=="0/0":
                anc=ref
                mut=alt
            elif gts[anc_sample_idx]=="1/1":
                anc=alt
                mut=ref
            else:
                continue
            
            tnc=reference[chrom][(pos-2):(pos+1)].seq.upper()
                
            is_CpG=(anc=="C" and mut=="T" and tnc[2]=="G") or (anc=="G" and mut=="A" and tnc[0]=="C")
                
            if is_CpG:
                CpG.write(line)
            else:
                notCpG.write(line)
        else:
            raise Exception("Not a comment line, but haven't seen header yet\n"+line)

    CpG.close()
    notCpG.close()

################################################################################

if __name__=="__main__":
    options=parse_options()
    main(options)


