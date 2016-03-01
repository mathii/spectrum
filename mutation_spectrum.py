'''
Created on Mar 5, 2015
From a vcf, compute the mutation spectrum, relative to one of the samples
of private mutations (i.e. singletons), in a trinucleotide context

@author: mathii
'''

from __future__ import division, print_function
from pyfaidx import Fasta
import sys, getopt, gzip

BASES=["A", "C", "G", "T"]

##########################################################################################################

def parse_options():
    """
    vcf: vcf input
    ref: 
    """
    options ={"vcf":None, "ref":None, "ref_sample":None, "out":"results", "count":1, "AA_INFO":False }

    try:
        opts, args = getopt.getopt(sys.argv[1:], "v:r:s:o:c:p:a",
                                    ["vcf", "ref", "ref_sample", "out", "count", "AA_INFO"])

    except Exception as err:
        print( str(err), file=sys.stderr)
        sys.exit()

    for o, a in opts:
        if o in ["-v","--vcf"]:             options["vcf"] = a
        elif o in ["-r","--ref"]:           options["ref"] = a
        elif o in ["-s","--ref_sample"]:    options["ref_sample"] = a
        elif o in ["-o","--out"]:           options["out"] = a
        elif o in ["-c","--count"]:         options["count"] = int(a) #allele count of variants to include
        elif o in ["-a", "AA_INFO"]:        options["AA_INFO"]=True

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

def make_dict(N_samples):
    '''
    Make the array to store the results - a dict with keys given by ('XYZ', 'W") which represents
    a trinucleotide mutation XYZ->XWZ
    '''
    
    d={}
    for a in BASES:
        for b in BASES:
            for c in BASES:
                for x in BASES:
                    if x!=b:
                        d[(a+b+c,x)]=[0]*N_samples
                         
    return d


##########################################################################################################

def print_results(results):
    out=open(options["out"], "w")
    out.write("Mutation\t"+"\t".join(results["samples"])+"\n")
    del(results["samples"])
    for k,v in results.items():
        out.write(k[0]+"."+k[1]+"\t"+"\t".join([str(x) for x in v])+"\n")


##########################################################################################################

def main(options):
    """
    Run the analysis. 
    """
    
    reference=Fasta(options["ref"])
    data=open2(options["vcf"])

    polarise_i=None
    results=None
    line_i=0
    skipped=0
    counted=0

    for line in data:
        line_i+=1
        if line.startswith("##"):
            continue
        elif line.startswith("#"): #Header line
            bits=line.split()
            N_samples=len(bits)-9
            results=make_dict(N_samples)
            if options["ref_sample"] and not options["AA_INFO"]:
                polarise_i=bits.index(options["ref_sample"])-9
            results["samples"]=bits[9:]
            
        else: #An actual line of data
#            if not line_i % 1000:
#                print("\rline: "+str(line_i))
            bits=line.split()
            chrom=bits[0]
            pos=int(bits[1])
            anc=ref=bits[3]
            mut=alt=bits[4]
            hetgts=["01", "10"]
            mutgt="11"
            
            if len(anc)>1 or len(ref)>1:
                continue

            if options["AA_INFO"]:
                info=bits[7]
                info_bits=info.split(";")
                AA=None
                for ib in info_bits:
                    if ib[:2]=="AA":
                        AA=ib[3]
                        break
                if AA==ref:
                    pass
                elif AA==alt:
                    anc=alt
                    mut=ref
                    mutgt="00"
                else:
                    skipped+=1
                    continue
            elif polarise_i:
                if bits[9+polarise_i]=="0/0":
                    pass
                elif bits[9+polarise_i]=="1/1":
                    anc=alt
                    mut=ref
                    mutgt="00"
                else:
                    skipped+=1
                    continue
                            
            tnc=reference[chrom][(pos-2):(pos+1)].seq
            if tnc[1]!=ref:
                skipped+=1
                continue
            
            # Check for the condition, exactly options["count"] of the allele in the dataset. 
            gts=[g[0]+g[2] for g in bits[9:]]
            het_count=sum([g in hetgts for g in gts])
            hom_count=sum([g==mutgt for g in gts])
            total_count=het_count+2*hom_count
            if total_count==options["count"]:
                counted+=1
                which_is_het= [i for i, x in enumerate(gts) if x in hetgts]
                which_is_hom= [i for i, x in enumerate(gts) if x==mutgt]
                key=(tnc[0]+anc+tnc[2], mut)
                if key in results:
                    for igt in which_is_het:
                        results[key][igt]+=1
                    for igt in which_is_hom:
                        results[key][igt]+=2
                else:
                    skipped+=1
                    
    print("Skipped "+str(skipped)+"/"+str(line_i))
    print("Counted "+str(counted)+"/"+str(line_i))
    print_results(results)
                     
##########################################################################################################


if __name__ == '__main__':
    options,args=parse_options()
    main(options)
