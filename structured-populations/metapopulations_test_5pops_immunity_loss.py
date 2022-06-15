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


npop = 5 # number of populations
completely_connected = False

if npop == 5:
    n_populations = ['A', 'B', 'C', 'D', 'E']
    n_hosts = [100, 100, 100, 100, 100]
    pops_dist = dict(zip(n_populations,n_hosts)) 
    infected_pops = ['A', 'B', 'C', 'D', 'E']
elif npop == 10:
    n_populations = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
    n_hosts = [100, 100, 100, 100, 100, 100, 100, 100, 100, 100]
    pops_dist = dict(zip(n_populations,n_hosts))
    infected_pops = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']


# Partially connected metapopulation
# For 5         For 10
# A   B           C     
#  \ /        B   |  D
#   C           \ | /
#  / \      A -- [c] -- E       
# D   E             

# [c] = clustered population of 5 populations

# Completely connected metapopulation
# A --- B
# |\   /|
# |  C  |
# |/   \|
# D --- E

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


num_hosts = sum(n_hosts)
inf_frac = 0.5 # base infected fraction
               # (this ultimately doesn't matter as long as it's high)


model = Model()

# Create Setup
model.newSetup( # Now, we'll define our new setup:
    'setup_normal', preset='host-host', # Use default host-host parameters.
    possible_alleles='AB',
    num_loci=len(wt_seq),
    fitnessHost=fitnessLandscape,
    mean_inoculum_host=1,
    contact_rate_host_host=0.02, #5.75e-3,
    recovery_rate_host=5e-3,
    mutate_in_host=2e-2,
    recombine_in_host=0,
    immunity_acquisition_rate_host=1e-2,
        # rate at which immunity is acquired within infected individuals
    immunity_loss_rate_host=1e-2,
        # rate at which immunity is lost within infected individuals
    immunityWeightsHost=myImmuneWeights,
        # immunity function that evaluates effect of immune memory on pathogen
        # genomes.
    )

# Create populations
for i in range(len(n_populations)):
    model.newPopulation(
        'pop' + n_populations[i],
        'setup_normal',
        num_hosts=pops_dist[n_populations[i]]
        )

# Create links
if completely_connected is False:
    # A - C
    model.linkPopulationsHostMigration('popA','popC',2e-3)
    model.linkPopulationsHostMigration('popC','popA',2e-3)
    # B - C
    model.linkPopulationsHostMigration('popB','popC',2e-3)
    model.linkPopulationsHostMigration('popC','popB',2e-3)
    # D - C
    model.linkPopulationsHostMigration('popD','popC',2e-3)
    model.linkPopulationsHostMigration('popC','popD',2e-3)
    # E - C
    model.linkPopulationsHostMigration('popE','popC',2e-3)
    model.linkPopulationsHostMigration('popC','popE',2e-3)

elif completely_connected is True:
    # A - B
    model.linkPopulationsHostMigration('popA','popB',2e-3)
    model.linkPopulationsHostMigration('popB','popA',2e-3)
    # A - D
    model.linkPopulationsHostMigration('popA','popD',2e-3)
    model.linkPopulationsHostMigration('popD','popA',2e-3)
    # A - C
    model.linkPopulationsHostMigration('popA','popC',2e-3)
    model.linkPopulationsHostMigration('popC','popA',2e-3)
    # B - C
    model.linkPopulationsHostMigration('popB','popC',2e-3)
    model.linkPopulationsHostMigration('popC','popB',2e-3)
    # B - E
    model.linkPopulationsHostMigration('popB','popE',2e-3)
    model.linkPopulationsHostMigration('popE','popB',2e-3)
    # D - C
    model.linkPopulationsHostMigration('popD','popC',2e-3)
    model.linkPopulationsHostMigration('popC','popD',2e-3)
    # D - E
    model.linkPopulationsHostMigration('popD','popE',2e-3)
    model.linkPopulationsHostMigration('popE','popD',2e-3)
    # E - C
    model.linkPopulationsHostMigration('popE','popC',2e-3)
    model.linkPopulationsHostMigration('popC','popE',2e-3)

# Add pathogens
for i in range(len(infected_pops)):
    model.addPathogensToHosts(
            'pop' + infected_pops[i],
            { wt_seq : int(pops_dist[n_populations[i]]*inf_frac) }
            # initial number of infected hosts
            # this ultimately doesn't matter as long as it's high, it'll stabilize
            # if it's too low, you might either get extinction or an expansion
            # event that could increase diversity (and evolution of resistance)
            # in high transmission
            )

# Run and save
output = model.run(0,100,time_sampling=0)
saveas = 1 if completely_connected is True else 0
savename = {1 : 'complete', 0 : 'partial'}
data = model.saveToDataFrame(
    'tests/5pops_sims/{}/metapopulations_migration_5pops.csv'.format(savename[saveas]))

# Plot
graph = model.populationsPlot( # Plot infected hosts per population over time.
    'tests/5pops_sims/{}/metapopulations_migration_5pops.png'.format(savename[saveas]), 
    data,
    num_top_populations=5, # plot all 5 populations
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
    data, track_specific_sequences=['8 A, 0 B','7 A, 1 B','6 A, 2 B','5 A, 3 B','4 A, 4 B','3 A, 5 B','2 A, 6 B','1 A, 7 B','0 A, 8 B'],
    num_top_sequences=-1,
    )

graph_composition = model.compositionPlot(
    'tests/5pops_sims/{}/pathogen_comp_metapopulations_migration_5pops.png'.format(savename[saveas]), 
    data,
    composition_dataframe=comp_dat,
    track_specific_sequences=['8 A, 0 B','7 A, 1 B','6 A, 2 B','5 A, 3 B','4 A, 4 B','3 A, 5 B','2 A, 6 B','1 A, 7 B','0 A, 8 B'],
    palette=CB_PALETTE_mod,
    stacked=False
    )

# Composition plot for each population:
for i, pop in enumerate(n_populations):
    graph_composition = model.compositionPlot(
    'tests/5pops_sims/{}/individual_pop/pathogen_comp_metapopulations_migration_pop{}.png'.format(savename[saveas], pop), 
    data,
    composition_dataframe=comp_dat,
    populations=['pop' + pop],
    track_specific_sequences=['8 A, 0 B','7 A, 1 B','6 A, 2 B','5 A, 3 B','4 A, 4 B','3 A, 5 B','2 A, 6 B','1 A, 7 B','0 A, 8 B'],
    palette=CB_PALETTE_mod,
    stacked=False
    )


# Gropup genomes by number of B alleles:
for genome in model.global_trackers['genomes_seen']:
    data['Immunity'] = data['Immunity'].str.replace(
        genome,
        str( len(genome)-genome.count('B') ) + ' A, ' + str(genome.count('B')) + ' B'
        )

# Plot genomes grouped by number of B alleles:
comp_dat = compositionDf(
    data, track_specific_sequences=['8 A, 0 B','7 A, 1 B','6 A, 2 B','5 A, 3 B','4 A, 4 B','3 A, 5 B','2 A, 6 B','1 A, 7 B','0 A, 8 B'],
    num_top_sequences=-1,
    )

graph_composition = model.compositionPlot(
    'tests/5pops_sims/{}/immunity_comp_metapopulations_migration_5pops.png'.format(savename[saveas]), 
    data,
    type_of_composition='Immunity', y_label='Genomes in immune memory',
    track_specific_sequences=['8 A, 0 B','7 A, 1 B','6 A, 2 B','5 A, 3 B','4 A, 4 B','3 A, 5 B','2 A, 6 B','1 A, 7 B','0 A, 8 B'],
    palette=CB_PALETTE_mod,
    stacked=False
    )

# Composition plot for each population:
for i, pop in enumerate(n_populations):
    graph_composition = model.compositionPlot(
    'tests/5pops_sims/{}/individual_pop/immunity_comp_metapopulations_migration_pop{}.png'.format(savename[saveas], pop), 
    data,
    populations=['pop' + pop],
    type_of_composition='Immunity', y_label='Genomes in immune memory',
    track_specific_sequences=['8 A, 0 B','7 A, 1 B','6 A, 2 B','5 A, 3 B','4 A, 4 B','3 A, 5 B','2 A, 6 B','1 A, 7 B','0 A, 8 B'],
    palette=CB_PALETTE_mod,
    stacked=False
    )