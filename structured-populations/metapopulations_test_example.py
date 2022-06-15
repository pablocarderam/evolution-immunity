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


model = Model()
model.newSetup('setup_normal', preset='vector-borne')
model.newSetup(
    'setup_cluster',
    contact_rate_host_vector = ( 2 *
        model.setups['setup_normal'].contact_rate_host_vector ),
    preset='vector-borne'
    ) # uses default parameters but doubles contact rate of the first setup

model.newPopulation('population_A','setup_normal', num_hosts=20, num_vectors=20)
model.newPopulation('population_B','setup_normal', num_hosts=20, num_vectors=20)
    # Create two populations that will be connected.
model.newPopulation(
    'isolated_population','setup_normal', num_hosts=20, num_vectors=20
    ) # A third population will remain isolated.

model.createInterconnectedPopulations(
    5,'clustered_population_','setup_cluster',
    host_migration_rate=2e-3, vector_migration_rate=0,
    host_host_contact_rate=0, vector_host_contact_rate=0,
    num_hosts=20, num_vectors=20
    )
    # Create a cluster of 5 populations connected to each other with a migration
    # rate of 2e-3 between each of them in both directions. Each population has
    # an numbered ID with the prefix "clustered_population_", has the parameters
    # defined in the "setup_cluster" setup, and has 20 hosts and vectors.
model.linkPopulationsHostMigration('population_A','clustered_population_4', )
    # We link population_A to one of the clustered populations with a one-way
    # migration rate of 2e-3.
model.linkPopulationsHostMigration('population_A','population_B',2e-3)
    # We link population_A to population_B with a one-way migration rate of
    # 2e-3.

model.addPathogensToHosts( 'population_A',{'AAAAAAAAAA':5} )
    # population_A starts with AAAAAAAAAA genotype pathogens.
model.addPathogensToHosts( 'population_B',{'GGGGGGGGGG':5} )
    # population_B starts with GGGGGGGGGG genotype pathogens.

output = model.run(0,100,time_sampling=0)
data = model.saveToDataFrame('tests/metapopulations_migration_example.csv')
graph = model.populationsPlot( # Plot infected hosts per population over time.
    'tests/metapopulations_migration_example.png', data,
    num_top_populations=8, # plot all 8 populations
    track_specific_populations=['isolated_population'],
        # Make sure to plot th isolated population totals if not in the top
        # infected populations.
    y_label='Infected hosts' # change y label
    )