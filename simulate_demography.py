# Simulate the demography of Africa, Eurasia and America
# And the frequency spectra of CpG and nonCpG mutations.

from __future__ import division,print_function
import msprime, math, gzip, random, pdb
import numpy as np

def out_of_africa():
    # First we set out the maximum likelihood values of the various parameters
    # given in Table 1.
    N_A = 12000
    N_B = 1000
    N_AF = 12000
    N_EU0 = 1000
    N_AS0 = 1000
    N_AM0 = 100
    # Times are provided in years, so we convert into generations.
    generation_time = 29
    T_AF = 70e3 / generation_time
    T_B = 60e3 / generation_time
    T_EU_AS = 40e3 / generation_time
    T_AS_AM = 20e3 / generation_time 
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU = 0.004
    r_AS = 0.0055
    r_AM = 0.01
    N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    N_AM = N_AM0 / math.exp(-r_AM * T_AS_AM)
    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB 3=AMerican
    # initially.
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=100, initial_size=N_AF),
        msprime.PopulationConfiguration(
            sample_size=100, initial_size=N_EU, growth_rate=r_EU),
        msprime.PopulationConfiguration(
            sample_size=100, initial_size=N_AS, growth_rate=r_AS),
        msprime.PopulationConfiguration(
            sample_size=40, initial_size=N_AM, growth_rate=r_AM)
    ]
    demographic_events = [
        # CEU and CHB merge into B with rate changes at T_EU_AS
        msprime.MassMigration(
            time=T_AS_AM, source=3, destination=2, proportion=1.0),
        msprime.MassMigration(
            time=T_EU_AS, source=2, destination=1, proportion=1.0),
        msprime.PopulationParametersChange(
            time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1),
        # Population B merges into YRI at T_B
        msprime.MassMigration(
            time=T_B, source=1, destination=0, proportion=1.0),
        # Size changes to N_A at T_AF
        msprime.PopulationParametersChange(
            time=T_AF, initial_size=N_A, population_id=0)
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    # dp = msprime.DemographyDebugger(
    #     Ne=N_A,
    #     population_configurations=population_configurations,
    #     demographic_events=demographic_events)
    # dp.print_history()
    return population_configurations, demographic_events

CPG_double_rate=0.25
SEQ_LEN=1e6
MU_SCALE=100
if __name__=="__main__":
    population_configurations, demographic_events =  out_of_africa()
    print("Simulating nonCpG")
    tree_sequence = msprime.simulate(
        population_configurations=population_configurations, 
        demographic_events=demographic_events,
        length=SEQ_LEN,
        mutation_rate=2e-8*MU_SCALE)

    print("Getting mutations nonCpG")
    # nonCpG_file=gzip.open("/Users/mathii/spectrum/sims/nonCpG_spectrum.gz", "w")
    shape = tree_sequence.get_num_mutations(), tree_sequence.get_sample_size()
    A = np.empty(shape, dtype="u1")
    for variant in tree_sequence.variants():
        A[variant.index] = variant.genotypes
    print("Printing nonCpG")
    # np.savetxt("/Users/mathii/spectrum/sims/nonCpG_spectrum.gz", A, delimiter="", fmt="%d")
    np.save("/Users/mathii/spectrum/sims/nonCpG_spectrum.npy", A)

    # for variant in tree_sequence.variants():
    #     # print(variant.index, variant.position, variant.genotypes, sep="\t")
    #     nonCpG_file.write("".join([str(x) for x in variant.genotypes])+"\n")

    print("Simulating CpG")
    tree_sequence = msprime.simulate(
        population_configurations=population_configurations,
        demographic_events=demographic_events,
        length=SEQ_LEN,
        mutation_rate=2e-8*MU_SCALE)

    print("Getting mutations CpG")
    # CpG_file=gzip.open("/Users/mathii/spectrum/sims/CpG_spectrum.gz", "w")
    shape = tree_sequence.get_num_mutations(), tree_sequence.get_sample_size()
    A = np.zeros(shape, dtype="u1")
    vs=tree_sequence.variants()
    for variant in vs:
        gt=variant.genotypes
        if random.random()<CPG_double_rate:
            try:
                v2=vs.next().genotypes
            except StopIteration:
                v2=gt
            A[variant.index]=np.array(np.logical_or(gt, v2), dtype=np.uint8)
    print("Printing CpG")
    # np.savetxt("/Users/mathii/spectrum/sims/CpG_spectrum.gz", A, delimiter="", fmt="%d")
    np.save("/Users/mathii/spectrum/sims/CpG_spectrum.npy", A)
