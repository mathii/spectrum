# Split a vcf accoriding to sites that are polymorphic in a specific population
#Default Africa, or not. 

from __future__ import division, print_function
import argparse, re, gzip, pdb

################################################################################

def parse_options():
    """
    argparse
    """
    parser=argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', type=str, default="")
    parser.add_argument('-m', '--map', type=str, default="")    
    parser.add_argument('-p', '--pop', type=str, default="Africa")
    return parser.parse_args()

################################################################################

def load_panel(panel_file):
    """
    indiviual->population map
    """
    panel={}
    pfile=open(panel_file, "r")
    for line in pfile:
        bits=line[:-1].split()
        panel[bits[0]]=bits[1]
    return panel

################################################################################

def main(options):
    
    panel=load_panel(options.map)
    
    is_gz=options.vcf[-3:]==".gz"
    
    basename=re.sub("\.vcf.+$", "", options.vcf)
    out_poly=basename+".poly_"+options.pop+".vcf"
    out_notpoly=basename+".notpoly_"+options.pop+".vcf"
    if is_gz:
        out_poly=out_poly+".gz"
        out_notpoly=out_notpoly+".gz"
        poly=gzip.open(out_poly, "w")
        notpoly=gzip.open(out_notpoly, "w")
        infile=gzip.open(options.vcf, "r")
    else:
        poly=open(out_poly, "w")
        notpoly=open(out_notpoly, "w")
        infile=open(options.vcf, "r")

    reading_data=False
    samples=[]
    for line in infile:        
        if not reading_data and line[:2]=="##":
            poly.write(line)
            notpoly.write(line)
        elif not reading_data and line[0]=="#":
            poly.write(line)
            notpoly.write(line)
            bits=line.split()
            samples=bits[9:]
            include_sample_idx=[i for i,x in enumerate(samples)if panel.get(x, "")==options.pop]
            reading_data=True
        elif reading_data:
            gt=[x[:3] for x in line.split()[9:]]
            pop_gt=[gt[i] for i in include_sample_idx]
            is_poly=not all([x==pop_gt[0] for x in pop_gt if x!="./."])
            if is_poly:
                poly.write(line)
            else:
                notpoly.write(line)
        else:
            raise Exception("Not a comment line, but haven't seen header yet\n"+line)

    poly.close()
    notpoly.close()

################################################################################

if __name__=="__main__":
    options=parse_options()
    main(options)

