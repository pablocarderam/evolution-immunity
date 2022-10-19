import sys
import os
import os.path
sys.path.append(
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt # plots
import seaborn as sns # pretty plots
import config as cf
from opqua.model import Model
from opqua.internal.data import compositionDf, compartmentDf
assert cf

path_save = 'tentative_figs/fig3_immunity_effect_over_evolution/immunity_effect_over_evolution_weights_dist_no_recom'

CB_PALETTE_mod = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0a2200",
                  "#0072B2", "#D55E00", "#f5d9ff", "#CC79A7", "#999999"]

num_hosts = 1000
wt_seq = 'AAAAA' # wild-type genome sequence
re_seq = 'BBBBB' # drug resistant mutant genome sequence
len_seq = 5
inf_frac = 0.5 # base infected fraction
               # (this ultimately doesn't matter as long as it's high)

recombination = True
recombination_rate = 2e-1 if recombination is True else 0

w1 = [0.4,0.25,0.15,0.1,0.1]  # decaying weights
w2 = [0.1,0.1,0.15,0.25,0.4]  # increasing weights
w3 = [0.1,0.2,0.4,0.2,0.1]    # bell weights
w4 = [0.4,0.05,0.1,0.05,0.4]  # U weights (ultra high)
w5 = [0.3,0.15,0.1,0.15,0.3]  # U weights (ultra soft)
pmw1 = [1]      # Perfect matching weight
pmw2 = [0.5]    # Perfect matching weight 

exps = {
    'decaying' : w1,
    'increasing' : w2,
    'bell' : w3,
    'U_high' : w4,
    'U_soft' : w5,
    'perfect_matching_all' : pmw1,
    'perfect_matching_half' : pmw2,
}

#### RUN MODEL ####

for key, elem in exps.items():    
    def myImmuneWeights(genome,immune_seq,weights=elem):
        if len(weights) == 1:
            return Model.perfectMatchImmunity(
                genome, immune_seq, weight=weights[0]
                )
            # Perfect matching: any mismatches between pathogen genome and immune
            # memory generate no immunity. Only perfect matches generate 100%
            # immunity.
        else:
            return Model.elementMatchImmunity(
                genome, immune_seq, weights
                )
            # Element matching: any mismatches  
    model = Model()
    model.newSetup( # Now, we'll define our new setup:
        'normal_setup', preset='host-host', # Use default host-host parameters.
        possible_alleles='AB',
            # Define "letters" in the "genome", or possible alleles for each locus.
            # Each locus can have different possible alleles if you define this
            # argument as a list of strings, but here, we take the simplest
            # approach.
        num_loci=len(wt_seq),
            # Define length of "genome", or total number of alleles.
        # recombinationHost= recombination_rate
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
        recombine_in_host = 0
        )

    model.newPopulation(
        'normal_pop','normal_setup', num_hosts=int( num_hosts )
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


    t_f = 15000 # final timepoint
    model.run(0,1000,time_sampling=100)

    if not os.path.isdir(os.path.join(path_save,key)):
        os.makedirs(os.path.join(path_save,key))

    data = model.saveToDataFrame(
        os.path.join(path_save,key,'weights_dist.csv')
        )

    graph_populations = model.populationsPlot( # Plot infected hosts per population over time.
        os.path.join(path_save,key,'populations_plot.png'), 
        data,
        num_top_populations=2, # plot all 2 populations
        y_label='Infected hosts' # change y label
        )
    
    # Gropup genomes by number of B alleles:
    for genome in model.global_trackers['genomes_seen']:
        data['Pathogens'] = data['Pathogens'].str.replace(
            genome,
            str( len(genome)-genome.count('B') ) + ' A, ' + str(genome.count('B')) + ' B'
            )

    # Plot genomes grouped by number of B alleles:
    comp_dat = compositionDf(
        data, track_specific_sequences=['5 A, 0 B','4 A, 1 B','3 A, 2 B','2 A, 3 B','1 A, 4 B','0 A, 5 B'],
        num_top_sequences=-1,
        )

    graph_composition = model.compositionPlot(
            # Create a plot to track pathogen genotypes across time.
        os.path.join(path_save,key,'pathogen_composition.png'),
        data,
        composition_dataframe=comp_dat,
        track_specific_sequences=['5 A, 0 B','4 A, 1 B','3 A, 2 B','2 A, 3 B','1 A, 4 B','0 A, 5 B'],
        type_of_composition='Pathogens',
        palette=CB_PALETTE_mod,
        # stacked=False
        )

    # Gropup genomes by number of B alleles:
    for genome in model.global_trackers['genomes_seen']:
        data['Immunity'] = data['Immunity'].str.replace(
            genome,
            str( len(genome)-genome.count('B') ) + ' A, ' + str(genome.count('B')) + ' B'
            )

    # Plot genomes grouped by number of B alleles:
    comp_dat = compositionDf(
        data, track_specific_sequences=['5 A, 0 B','4 A, 1 B','3 A, 2 B','2 A, 3 B','1 A, 4 B','0 A, 5 B'],
        num_top_sequences=-1,
        )

    graph_composition = model.compositionPlot(
            # Create a plot of genotypes in the hosts' immune memories across time.
        os.path.join(path_save,key,'immunity_composition.png'), 
        data,
        # composition_dataframe=comp_dat,
        type_of_composition='Immunity', y_label='Genomes in immune memory',
        track_specific_sequences=['5 A, 0 B','4 A, 1 B','3 A, 2 B','2 A, 3 B','1 A, 4 B','0 A, 5 B'],
        palette=CB_PALETTE_mod,
        # stacked=False
        )

    graph_compartments = model.compartmentPlot(
        os.path.join(path_save,key,'reassortment_compartments.png'), 
        data
        )
        # Also generate a normal compartment plot. Notice the total number of
        # infections in the composition plot can exceed the number of infected hosts
        # in the compartment plot. This happens because a single host infected by
        # multiple genotypes is counted twice in the former, but not the latter.

    # Plot all genomes separately:
    graph_compartments = model.compartmentPlot(
        os.path.join(path_save,key,'weights_dist_compartments.png'), data,
        # track_specific_sequences=['5 A, 0 B','4 A, 1 B','3 A, 2 B','2 A, 3 B','1 A, 4 B','0 A, 5 B'],
        )

    graph_composition = model.compositionPlot(
        os.path.join(path_save,key,'weights_dist_composition.png'), data,
        # track_specific_sequences=['5 A, 0 B','4 A, 1 B','3 A, 2 B','2 A, 3 B','1 A, 4 B','0 A, 5 B'],
        count_individuals_based_on_model=model
        )

    graph_clustermap = model.clustermap(
        os.path.join(path_save,key,'weights_dist_clustermap.png'), data,
        save_data_to_file='weights_dist_pairwise_distances.csv',
        num_top_sequences=15,
        )

    # Plot weights distribution:
    df = pd.DataFrame()
    df['pos'] = [0,1,2,3,4]
    df['weight'] = elem
    fig, ax = plt.subplots(1, 1, figsize=(5,5))
    sns.barplot(ax=ax, data=df,x='pos', y='weight', label=key)
    #eliminate x axis
    ax.set_xlabel('')
    ax.set_ylim(0,1.01)
    # ax.set_xticks([])
    ax.set_ylabel('Weight')
    fig.savefig(os.path.join(path_save,key,'weights_dist.png'))