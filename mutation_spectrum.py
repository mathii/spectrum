'''
Created on Mar 5, 2015
From a vcf, compute the mutation spectrum, relative to one of the samples
of private mutations (i.e. singletons), in a trinucleotide context

@author: mathii
'''

from __future__ import division, print_function
from pyfaidx import Fasta
import sys, getopt, gzip
from collections import defaultdict

BASES=["A", "C", "G", "T"]

##########################################################################################################

def parse_options():
    """
    vcf: vcf input
    ref: 
    """
    options ={"vcf":None, "ref":None, "ref_sample":None, "out":"results", "count":None,
              "AA_INFO":False, "mpf":None, "filter_file":None, "filter_values":(),
              "filter_list":None, "pos_out":None, "private_panel":None, "hom_weight":2,
              "max_missing":0, "n_context":3}

    try:
        opts, args = getopt.getopt(sys.argv[1:], "v:r:s:o:m:c:p:n:f:l:i:h:x:at:",
        ["vcf", "ref", "ref_sample", "out", "mpf", "count", "pos_out", "private_panel", "filter_file", "filter_value", "filter_list", "hom_weight", "max_missing", "AA_INFO", "n_context"])

    except Exception as err:
        print( str(err), file=sys.stderr)
        sys.exit()

    for o, a in opts:
        if o in ["-v","--vcf"]:             options["vcf"] = a
        elif o in ["-r","--ref"]:           options["ref"] = a
        elif o in ["-s","--ref_sample"]:    options["ref_sample"] = a
        elif o in ["-o","--out"]:           options["out"] = a
        elif o in ["-h","--hom_weight"]:    options["hom_weight"] = float(a)
        elif o in ["-x","--max_missing"]:   options["max_missing"] = int(a)
        elif o in ["-m","--mpf"]:           options["mpf"] = a #Output mutation position format
        elif o in ["-p","--pos_out"]:       options["pos_out"] = a #Output position format with all variants and contex.
        elif o in ["-n","--private_panel"]: options["private_panel"] = a   #Only output private mutations
        elif o in ["-c","--count"]:         options["count"] = int(a) #allele count of variants to include
        elif o in ["-a", "AA_INFO"]:        options["AA_INFO"]=True
        elif o in ["-f","--filter_file"]:   options["filter_file"] = a #Filter values according to this fasta file
        elif o in ["-l","--filter_value"]:  options["filter_values"] = set(a.split(",")) #Include sites that match these values in the fasta
        elif o in ["-i","--filter_list"]:   options["filter_list"] = a #A file containing a specific list of sites to include chr/pos
        elif o in ["-t","--n_context"]:     options["n_context"] = int(a) #Nucleotide context
    print( "found options:", file=sys.stderr)
    print( options, file=sys.stderr)

    if not options["n_context"] % 2 or options["n_context"]<1:
        raise Exception("nucleotide context must be a positive odd number")

    if not options["count"]:
        print("***** Using variants of all counts ******")
    else:
        print("Using variants of derived allele count "+str(options["count"]))
            
    
    return options, args

##########################################################################################################

def open2(f, mode="r"):
    if f.endswith(".gz"):
        return gzip.open(f, mode)
    else:
        return open(f, mode)

##########################################################################################################

def make_dict(N_samples, options):
    '''
    Make the array to store the results - a dict with keys given by ('XYZ', 'W") which represents
    a trinucleotide mutation XYZ->XWZ
    '''
    
    d={}
    keys=[""]
    for i in range(options["n_context"]):
        new_keys=[]
        for a in BASES:
            for k in keys:
                new_keys.append(a+k)
        keys=new_keys
        
    midpoint=options["n_context"]//2
    for k in keys:
        for a in BASES:
            if a!=k[midpoint]:
                d[(k,a)]=[0]*N_samples
                
    return d


##########################################################################################################

def read_filter_list(filter_list_file):
    """
    Two columns, chrom, pos of sites that should be included. 
    """
    filter_list=defaultdict(set)
    
    list_file=open2(filter_list_file, "r")
    for line in list_file:
        if line[0]=="#":
            continue
        bits=line.split()
        filter_list[bits[0]].add(int(bits[1]))
    
    return filter_list

##########################################################################################################

def print_results(results, options):
    out=open(options["out"], "w")
    out.write("Mutation\t"+"\t".join(results["samples"])+"\n")
    del(results["samples"])
    for k,v in results.items():
        out.write(k[0]+"."+k[1]+"\t"+"\t".join([str(x) for x in v])+"\n")


##########################################################################################################

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
        
##########################################################################################################

def main(options):
    """
    Run the analysis. 
    """
    
    reference=Fasta(options["ref"])
    data=open2(options["vcf"])

    filter=None
    if options["filter_file"]:
        filter=Fasta(options["filter_file"])

    filter_list=None
    if options["filter_list"]:
        filter_list=read_filter_list(options["filter_list"])

    if options["mpf"]:
        mpf_out=open(options["mpf"], "w")

    if options["pos_out"]:
        pos_out=open(options["pos_out"], "w")

    private_panel=None
    if options["private_panel"]:
        private_panel=load_panel(options["private_panel"])
        
    polarise_i=None
    results=None
    line_i=0
    skipped=0
    missing=0
    counted=0

    half_n_context=options["n_context"]//2

    for line in data:
        line_i+=1
        if line.startswith("##"):
            continue
        elif line.startswith("#"): #Header line
            bits=line.split()
            N_samples=len(bits)-9
            results=make_dict(N_samples, options)
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
            
            if len(alt)>1 or len(ref)>1: #Only include bialleleic SNPs
                continue

            if filter and filter[chrom][pos-1].seq not in options["filter_values"]:
                continue 

            if filter_list and pos not in filter_list[chrom]:
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
            elif options["ref_sample"]:
                if bits[9+polarise_i]=="0/0":
                    pass
                elif bits[9+polarise_i]=="1/1":
                    anc=alt
                    mut=ref
                    mutgt="00"
                else:
                    skipped+=1
                    continue
                            
            tnc=reference[chrom][(pos-1-half_n_context):(pos+half_n_context)].seq.upper()
            if tnc[half_n_context]!=ref:
                skipped+=1
                continue
            
            # Check for the condition, exactly options["count"] of the allele in the dataset.

            gtbits=bits[9:]
            if options["ref_sample"]:
                gtbits.pop(polarise_i)
            gts=[g[0]+g[2] for g in gtbits]   #Genotypes excludeing the polarising sample
            gts_with_pol=[g[0]+g[2] for g in bits[9:]]   #Genotypes with the polarising sample
            het_count=sum([g in hetgts for g in gts])
            hom_count=sum([g==mutgt for g in gts])
            miss_count=sum([g=="./." for g in gts])
            total_count=het_count+2*hom_count

            include_this_mutation=True
            if options["count"] and total_count!=options["count"]:
                include_this_mutation=False

            if options["max_missing"] and miss_count>options["max_missing"]:
                missing+=1
                include_this_mutation=False

            which_is_het= [i for i, x in enumerate(gts_with_pol) if x in hetgts]
            which_is_hom= [i for i, x in enumerate(gts_with_pol) if x==mutgt]
                
            if include_this_mutation and private_panel:
                pops=[private_panel[results["samples"][x]] for x in which_is_het+which_is_hom if results["samples"][x] != options["ref_sample"] ]
                if not len(pops) or not all([p==pops[0] for p in pops]):
                    include_this_mutation=False

            if include_this_mutation:
                counted+=1
                tac=list(tnc)
                tac[half_n_context]=anc
                tac="".join(tac)
                key=(tac, mut)
                if key in results:
                    for igt in which_is_het:
                        results[key][igt]+=1
                    for igt in which_is_hom:
                        results[key][igt]+=options["hom_weight"]
                else:
                    skipped+=1
                    
                if options["pos_out"]:
                    pos_out.write(chrom+"\t"+str(pos)+"\t"+tnc+"."+mut+"\n")
                    
                if options["mpf"]:
                    for hi in which_is_het:
                        if hi!=polarise_i:
                            mpf_out.write(results["samples"][hi]+"\tchr"+chrom+"\t"+str(pos)+"\t"+anc+"\t"+mut+"\n")
                    for hi in which_is_hom:
                        if hi!=polarise_i:
                            mpf_out.write(results["samples"][hi]+"\tchr"+chrom+"\t"+str(pos)+"\t"+anc+"\t"+mut+"\n")
                            mpf_out.write(results["samples"][hi]+"\tchr"+chrom+"\t"+str(pos)+"\t"+anc+"\t"+mut+"\n")                        
                    
    print("Skipped "+str(skipped)+"/"+str(line_i))
    print("Missing "+str(missing)+"/"+str(line_i))
    print("Counted "+str(counted)+"/"+str(line_i))
    print_results(results, options)
                     
##########################################################################################################


if __name__ == '__main__':
    options,args=parse_options()
    main(options)
