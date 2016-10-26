# Read a vcf output heterozygosity per-sample for Ts and Tv

from __future__ import division, print_function
import gzip, fileinput

HETS=["0/1", "1/0", "0|1", "1|0"]
MISS=["./.", ".|."]
TS=["AG", "CT"]
TV=["AC", "AT", "CG", "GT"] 

samples=[]
results={}
reading_data=False
for line in fileinput.input():
    if reading_data:
        bits=line[:-1].split()
        alleles="".join(sorted(bits[3:5]))
        for i,gt in enumerate(bits[9:]):
            if gt not in MISS:
                results[samples[i]]["Total"]+=1
        if alleles in TS:
            for i,gt in enumerate(bits[9:]):
                if gt in HETS:
                    results[samples[i]]["Ts"]+=1
        elif alleles in TV:
            for i,gt in enumerate(bits[9:]):
                if gt in HETS:
                    results[samples[i]]["Tv"]+=1
        else:
            for i,gt in enumerate(bits[9:]):
                if gt in HETS:
                    results[samples[i]]["Other"]+=1

    elif line[:6]=="#CHROM":
        samples=line[:-1].split()[9:]
        for s in samples:
            results[s]={"Ts":0, "Tv":0, "Other":0, "Total":0}
        reading_data=True
        
for s in samples:
    print("\t".join([s]+[str(results[s][x]) for x in ["Ts", "Tv", "Other", "Total"]]))
