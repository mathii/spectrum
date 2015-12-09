'''
Created on Mar 5, 2015
From a vcf, compute the mutation spectrum, relative to one of the samples
of CC->TT dinucleotide mutations
@author: mathii
'''

from __future__ import division, print_function
import sys, getopt, gzip, pdb

##########################################################################################################

def parse_options():
    """
    vcf: vcf input
    ref: 
    """
    options ={"vcf":None, "ref_sample":None, "count":1,  }

    try:
        opts, args = getopt.getopt(sys.argv[1:], "v:s:o:c:p:",
                                    ["vcf", "ref_sample", "out", "count"])

    except Exception as err:
        print( str(err), file=sys.stderr)
        sys.exit()

    for o, a in opts:
        if o in ["-v","--vcf"]:             options["vcf"] = a
        elif o in ["-s","--ref_sample"]:    options["ref_sample"] = a
        elif o in ["-c","--count"]:         options["count"] = int(a) #allele count of variants to include
    
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

def main(options):
    """
    Run the analysis. 
    """
    
    data=open2(options["vcf"])

    polarise_i=None
    results=None
    line_i=0
    skipped=0
    
    
    for line in data:
        line_i+=1
        if line.startswith("##"):
            continue
        elif line.startswith("#"): #Header line
            bits=line.split()
            N_samples=len(bits)-9
            if options["ref_sample"]:
                polarise_i=bits.index(options["ref_sample"])-9
            results_samples=bits[9:]
            results=[0]*N_samples
        else: #An actual line of data
            if last_line.startswith("#"):
                last_line=line
                continue

            bits=line.split()
            last_bits=last_line.split()
            
            pos=int(bits[1])
            last_pos=int(last_bits[1])

            if pos!=last_pos+1:
                last_line=line
                continue
                
            anc=ref=bits[3]
            mut=alt=bits[4]

            last_anc=last_ref=last_bits[3]
            last_mut=last_alt=last_bits[4]

            hetgts=["0/1", "1/0"]
            mutgt=last_mutgt="1/1"
            
            if polarise_i:
                if bits[9+polarise_i]=="0/0":
                    pass
                elif bits[9+polarise_i]=="1/1":
                    anc=alt
                    mut=ref
                    mutgt="0/0"
                else:
                    skipped+=1
                    last_line=line
                    continue

                if last_bits[9+polarise_i]=="0/0":
                    pass
                elif bits[9+polarise_i]=="1/1":
                    last_anc=alt
                    last_mut=ref
                    last_mutgt="0/0"
                else:
                    skipped+=1
                    last_line=line
                    continue

            if last_anc!="C" or last_mut!="T" or anc!="C" or mut!="T":
                    last_line=line
                    continue

            # Check for the condition, exactly [count alleles] and everyone else ancestral
            gts=bits[9:]
            last_gts=bits[9:]
            hets=[g in hetgts for g in gts]
            homs=[g==mutgt for g in gts]
            last_hets=[g in hetgts for g in last_gts]
            last_homs=[g==mutgt for g in last_gts]

            match=all([x==y for x,y in zip(hets, last_hets)]) and all([x==y for x,y in zip(homs, last_homs)])

            if not match:
                last_line=line
                continue
            
            total_count=sum(hets)+sum(homs)
            if total_count==options["count"]:
                which_is_het= [i for i, x in enumerate(gts) if x in hetgts]
                which_is_hom= [i for i, x in enumerate(gts) if x==mutgt]
                for igt in which_is_het:
                    results[igt]+=1
                for igt in which_is_hom:
                    results[igt]+=2


        last_line=line

    print("Skipped "+str(skipped)+"/"+str(line_i), file=sys.stdout)
    for i,x in enumerate(results):
        print(results_samples[i]+"\t"+str(x))
    
##########################################################################################################


if __name__ == '__main__':
    options,args=parse_options()
    main(options)
