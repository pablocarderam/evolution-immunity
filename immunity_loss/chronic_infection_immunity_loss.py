import numpy as np
import pandas as pd
import matplotlib.pyplot as plt # plots
import seaborn as sns # pretty plots

from opqua.model import Model
from opqua.internal.data import compositionDf, compartmentDf

CB_PALETTE_mod = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0a2200",
                  "#0072B2", "#D55E00", "#f5d9ff", "#CC79A7", "#999999"]
    # www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
    # http://jfly.iam.u-tokyo.ac.jp/color/


wt_seq = 'AAA'

def myImmuneWeights(genome,immune_seq):
    return Model.perfectMatchImmunity(
        genome, immune_seq, weight=1
        )
        # Perfect matching: any mismatches between pathogen genome and immune
        # memory generate no immunity. Only perfect matches generate 100%
        # immunity.

num_hosts = 1000
chronic_frac = 0.10
inf_frac = 0.5 # base infected fraction
               # (this ultimately doesn't matter as long as it's high)

model = Model()
model.newSetup( # Now, we'll define our new setup:
    'normal_setup', preset='host-host', # Use default host-host parameters.
    possible_alleles='AB',
        # Define "letters" in the "genome", or possible alleles for each locus.
        # Each locus can have different possible alleles if you define this
        # argument as a list of strings, but here, we take the simplest
        # approach.
    num_loci=3,
        # Define length of "genome", or total number of alleles.
    mutate_in_host=1e-2,
        # Modify de novo mutation rate of pathogens when in host to get some
        # evolution!
    immunity_acquisition_rate_host=1e-2,
        # rate at which immunity is acquired within infected individuals
    immunity_loss_rate_host=0,
        # rate at which immunity is lost within infected individuals
    immunityWeightsHost=myImmuneWeights,
        # immunity function that evaluates effect of immune memory on pathogen
        # genomes.
    contact_rate_host_host=2e-2,
    recovery_rate_host=5e-3,
    )

model.newSetup( # Now, we'll define our new setup:
    'chronic_setup', preset='host-host', # Use default host-host parameters.
    possible_alleles='AB',
        # Define "letters" in the "genome", or possible alleles for each locus.
        # Each locus can have different possible alleles if you define this
        # argument as a list of strings, but here, we take the simplest
        # approach.
     num_loci=3,
        # Define length of "genome", or total number of alleles.
    mutate_in_host=1e-2,
        # Modify de novo mutation rate of pathogens when in host to get some
        # evolution!
    immunity_acquisition_rate_host=1e-2,
        # rate at which immunity is acquired within infected individuals
    immunity_loss_rate_host=0,
        # rate at which immunity is lost within infected individuals
    immunityWeightsHost=myImmuneWeights,
        # immunity function that evaluates effect of immune memory on pathogen
        # genomes.
    contact_rate_host_host=2e-2,
    recovery_rate_host=5e-2,
    )

model.newPopulation(
    'normal_pop','normal_setup', num_hosts=int( num_hosts * (1-chronic_frac) )
    )
model.newPopulation(
    'chronic_pop','normal_setup', num_hosts=int( num_hosts * chronic_frac )
    )

model.addPathogensToHosts(
    'normal_pop',
    { wt_seq : int(num_hosts*inf_frac) }
    # initial number of infected hosts
    # this ultimately doesn't matter as long as it's high, it'll stabilize
    # if it's too low, you might either get extinction or an expansion
    # event that could increase diversity (and evolution of resistance)
    # in high transmission
    )

model.linkPopulationsHostHostContact('normal_pop','chronic_pop',2e-2)
model.linkPopulationsHostHostContact('chronic_pop','normal_pop',2e-2)

t_f = 15000 # final timepoint
model.run(0,1000,time_sampling=0)
data = model.saveToDataFrame(
    'tests/immunity_loss_v2_0.1_chronic_no_imm_loss/chronic_infection_example.csv'
    )

graph_populations = model.populationsPlot( # Plot infected hosts per population over time.
    'tests/immunity_loss_v2_0.1_chronic_no_imm_loss/populations_plot.png', 
    data,
    num_top_populations=2, # plot all 2 populations
    y_label='Infected hosts' # change y label
    )

graph_composition = model.compositionPlot(
        # Create a plot to track pathogen genotypes across time.
    'tests/immunity_loss_v2_0.1_chronic_no_imm_loss/pathogen_composition.png',
    data,
    type_of_composition='Pathogens',
    num_top_sequences=8,
        # Track the 8 most represented genomes overall (only 8 possible).
    )

graph_composition = model.compositionPlot(
        # Create a plot of genotypes in the hosts' immune memories across time.
    'tests/immunity_loss_v2_0.1_chronic_no_imm_loss/immunity_composition.png', 
    data,
    type_of_composition='Immunity', y_label='Genomes in immune memory',
    num_top_sequences=8,
        # Track the 8 most represented genomes overall (only 8 possible).
    )

graph_compartments = model.compartmentPlot(
    'tests/immunity_loss_v2_0.1_chronic_no_imm_loss/reassortment_compartments.png', 
    data
    )
    # Also generate a normal compartment plot. Notice the total number of
    # infections in the composition plot can exceed the number of infected hosts
    # in the compartment plot. This happens because a single host infected by
    # multiple genotypes is counted twice in the former, but not the latter.

graph_composition_A = model.compositionPlot(
        # Create a plot to track pathogen genotypes across time.
    'tests/immunity_loss_v2_0.1_chronic_no_imm_loss/pathogen_composition_fitness_valley_NORMAL_POP.png',
    data,
    populations=['normal_pop'],
    type_of_composition='Pathogens',
    num_top_sequences=8,
        # Track the 8 most represented genomes overall (only 8 possible).
    )

graph_composition_A = model.compositionPlot(
        # Create a plot of genotypes in the hosts' immune memories across time.
    'tests/immunity_loss_v2_0.1_chronic_no_imm_loss/immunity_composition_fitness_valley_NORMAL_POP.png', 
    data,
    populations=['normal_pop'],
    type_of_composition='Immunity', y_label='Genomes in immune memory',
    num_top_sequences=8,
        # Track the 8 most represented genomes overall (only 8 possible).
    )

graph_composition_B = model.compositionPlot(
        # Create a plot to track pathogen genotypes across time.
    'tests/immunity_loss_v2_0.1_chronic_no_imm_loss/pathogen_composition_fitness_valley_CHRONIC_POP.png',
    data,
    populations=['chronic_pop'],
    type_of_composition='Pathogens',
    num_top_sequences=8,
        # Track the 8 most represented genomes overall (only 8 possible).
    )

graph_composition_B = model.compositionPlot(
        # Create a plot of genotypes in the hosts' immune memories across time.
    'tests/immunity_loss_v2_0.1_chronic_no_imm_loss/immunity_composition_fitness_valle y_CHRONIC_POP.png', 
    data,
    populations=['chronic_pop'],
    type_of_composition='Immunity', y_label='Genomes in immune memory',
    num_top_sequences=8,
        # Track the 8 most represented genomes overall (only 8 possible).
    )

# Plot all genomes separately:
graph_compartments = model.compartmentPlot(
    'tests/immunity_loss_v2_0.1_chronic_no_imm_loss/chronic_infection_example_compartments.png', data
    )

graph_composition = model.compositionPlot(
    'tests/immunity_loss_v2_0.1_chronic_no_imm_loss/chronic_infection_example_composition_ind.png', data,
    num_top_sequences=6,
    count_individuals_based_on_model=model
    )

graph_clustermap = model.clustermap(
    'tests/immunity_loss_v2_0.1_chronic_no_imm_loss/chronic_infection_example_clustermap.png', data,
    save_data_to_file='chronic_infection_example_pairwise_distances.csv',
    num_top_sequences=15,
    )