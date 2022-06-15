 
'''
Cardenas, Torres, & Santos-Vega, 2022
Coded by github.com/pablocarderam

'''

import textdistance as td

import numpy as np
import matplotlib.pyplot as plt # plots
import seaborn as sns # pretty plots

from opqua.model import Model
from opqua.internal.data import compositionDf, compartmentDf

CB_PALETTE_mod = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0a2200",
                  "#0072B2", "#D55E00", "#f5d9ff", "#CC79A7", "#999999"]
    # www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
    # http://jfly.iam.u-tokyo.ac.jp/color/

wt_seq = 'AAAAAAAA' # wild-type genome sequence
re_seq = 'BBBBBBBB' # drug resistant mutant genome sequence

# Defines intra-host competitive fitness of pathogens
def fitnessLandscape(genome):
    num_b = genome.count('B')
    if num_b < 3:
        f = 1 / ( 2 * ( (num_b+1)**5 ) )
    else:
        f = 1 / ( ( len(genome)-num_b+1 )**3 )

    return f

# Define custom immunity function for the host:
# Immunity functions must take in 2 arguments and return a 0-1 immunity
# coefficient value. Here, we take advantage of one of the preset functions,
# but you can define it any way you want!
def myImmuneWeights(genome,immune_seq):
    return Model.perfectMatchImmunity(
        genome, immune_seq, weight=1
        )
        # Perfect matching: any mismatches between pathogen genome and immune
        # memory generate no immunity. Only perfect matches generate 100%
        # immunity.

# Modifies recovery rate of resistant pathogens, e.g. if hosts are seeking
# treatment by a drug when they get sick and the resistant pathogen endures for
# longer.
# Included in example as a way of increasing mutant persistance, makes effect
# easier to visualize.
def drugCheck(genome,re_seq='BBBBBBBB'):
    return 1 if genome!=re_seq else 0.75

num_hosts = 1000
chronic_frac = 0.05
inf_frac = 0.5 # base infected fraction
               # (this ultimately doesn't matter as long as it's high)

t_d = 10000 # time at which drug is applied

model = Model()
model.newSetup( # Now, we'll define our new setup:
    'normal_setup', preset='host-host', # Use default host-host parameters.
    possible_alleles='AB',
    num_loci=len(wt_seq),
    fitnessHost=fitnessLandscape,
    recoveryHost=drugCheck,
    mean_inoculum_host=1,
    contact_rate_host_host=0.02, #5.75e-3,
    recovery_rate_host=5e-3,
    mutate_in_host=2e-2,
    recombine_in_host=0, 
    immunity_acquisition_rate_host=1e-2,
        # rate at which immunity is acquired within infected individuals
    immunity_loss_rate_host=0,
        # rate at which immunity is lost within infected individuals
    immunityWeightsHost=myImmuneWeights,
        # immunity function that evaluates effect of immune memory on pathogen
        # genomes.
    )

model.newSetup( # Now, we'll define our new setup:
    'chronic_setup', preset='host-host', # Use default host-host parameters.
    possible_alleles='AB',
    num_loci=len(wt_seq),
    fitnessHost=fitnessLandscape,
    recoveryHost=drugCheck,
    mean_inoculum_host=1,
    contact_rate_host_host=0.02, #5.75e-3,
    recovery_rate_host=5e-2,
    mutate_in_host=2e-2,
    recombine_in_host=0,
    immunity_acquisition_rate_host=1e-2,
        # rate at which immunity is acquired within infected individuals
    immunity_loss_rate_host=0,
        # rate at which immunity is lost within infected individuals
    immunityWeightsHost=myImmuneWeights,
        # immunity function that evaluates effect of immune memory on pathogen
        # genomes.
    )

model.newPopulation(
    'normal_pop','normal_setup', num_hosts=int( num_hosts * (1-chronic_frac) )
    )
model.newPopulation(
    'chronic_pop','normal_setup', num_hosts=int( num_hosts * chronic_frac )
    )
# model.linkPopulationsHostHostContact('normal_pop','chronic_pop',5.75e-3)

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

model.newIntervention( t_d, 'treatHosts', [ 'normal_pop', 1, [re_seq] ] )
    # add drug, kills everything except resistant mutants

model.newIntervention( t_d, 'treatHosts', [ 'chronic_pop', 1, [re_seq] ] )
    # add drug, kills everything except resistant mutants

t_f = t_d+5000 # final timepoint
model.run(0,t_f,time_sampling=999, host_sampling=0, vector_sampling=0)
data = model.saveToDataFrame(
    'tests/intervential_all_pop/chronic_infection_example.csv'
    )

# Plot all genomes separately:
graph_compartments = model.compartmentPlot(
    'tests/intervential_all_pop/chronic_infection_example_compartments.png', data
    )

graph_composition = model.compositionPlot(
    'tests/intervential_all_pop/chronic_infection_example_composition_ind.png', data,
    num_top_sequences=6,
    count_individuals_based_on_model=model
    )

graph_clustermap = model.clustermap(
    'tests/intervential_all_pop/chronic_infection_example_clustermap.png', data,
    save_data_to_file='chronic_infection_example_pairwise_distances.csv',
    num_top_sequences=15,
    )

# Gropup genomes by number of B alleles:
for genome in model.global_trackers['genomes_seen']:
    data['Pathogens'] = data['Pathogens'].str.replace(
        genome,
        str( len(genome)-genome.count('B') ) + ' A, ' + str(genome.count('B')) + ' B'
        )

# Plot genomes grouped by number of B alleles:
comp_dat = compositionDf(
    data, track_specific_sequences=['8 A, 0 B','7 A, 1 B','6 A, 2 B','5 A, 3 B','4 A, 4 B','3 A, 5 B','2 A, 6 B','1 A, 7 B','0 A, 8 B'],
    num_top_sequences=-1,
    )

graph_composition = model.compositionPlot(
    'tests/intervential_all_pop/chronic_infection_example_composition.png', data,
    composition_dataframe=comp_dat,
    track_specific_sequences=['8 A, 0 B','7 A, 1 B','6 A, 2 B','5 A, 3 B','4 A, 4 B','3 A, 5 B','2 A, 6 B','1 A, 7 B','0 A, 8 B'],
    palette=CB_PALETTE_mod
    )

graph_pop = model.populationsPlot( # Plot infected hosts per population over time.
    'tests/intervential_all_pop/chronic_infection_example_populations.png', data,
    y_label='Infected hosts' # change y label
    )
