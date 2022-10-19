# Tests using genome of length 3

## Init paramenters:

    num_hosts = 1000
    wt_seq = 'AAA' # wild-type genome sequence
    re_seq = 'BBB' # drug resistant mutant genome sequence
    len_seq = 3
    inf_frac = 0.5 # base infected fraction
                # (this ultimately doesn't matter as long as it's high)

## Weights distributions:

    w1 = [0.5,0.3,0.2]  # decaying weights
    w2 = [0.2,0.3,0.5]  # increasing weights
    w3 = [0.25,0.5,0.25]    # bell weights
    w4 = [0.4,0.2,0.4]  # U weights 
    pmw1 = [1]      # Perfect matching weight

<table><tr>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/decaying/weights_dist.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/increasing/weights_dist.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/bell/weights_dist.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/U/weights_dist.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/weights_dist_perfect_matching.png" alt="Drawing" style="width: 250px;"/> </td>
</tr></table>

##### *This is the same order for the figures below*

## First we'll try with **recombination turned off**

#### Immune memory
<table><tr>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/decaying/immunity_composition.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/increasing/immunity_composition.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/bell/immunity_composition.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/U/immunity_composition.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/perfect_matching_all/immunity_composition.png" alt="Drawing" style="width: 250px;"/> </td>
</tr></table>

#### Pathogen composition
<table><tr>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/decaying/pathogen_composition.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/increasing/pathogen_composition.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/bell/pathogen_composition.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/U/pathogen_composition.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/perfect_matching_all/pathogen_composition.png" alt="Drawing" style="width: 250px;"/> </td>
</tr></table>

#### Composition
<table><tr>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/decaying/weights_dist_composition.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/increasing/weights_dist_composition.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/bell/weights_dist_composition.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/U/weights_dist_composition.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/perfect_matching_all/weights_dist_composition.png" alt="Drawing" style="width: 250px;"/> </td>
</tr></table>

#### Compartments
<table><tr>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/decaying/reassortment_compartments.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/increasing/reassortment_compartments.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/bell/reassortment_compartments.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/U/reassortment_compartments.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/perfect_matching_all/reassortment_compartments.png" alt="Drawing" style="width: 250px;"/> </td>
</tr></table>

#### Infections
<table><tr>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/decaying/populations_plot.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/increasing/populations_plot.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/bell/populations_plot.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/U/populations_plot.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/perfect_matching_all/populations_plot.png" alt="Drawing" style="width: 250px;"/> </td>
</tr></table>

#### Cluster map
<table><tr>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/decaying/weights_dist_clustermap.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/increasing/weights_dist_clustermap.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/bell/weights_dist_clustermap.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/U/weights_dist_clustermap.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_off/perfect_matching_all/weights_dist_clustermap.png" alt="Drawing" style="width: 250px;"/> </td>
</tr></table>

## Now with **recombination turned ON**
****
#### Immune memory
<table><tr>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/decaying/immunity_composition.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/increasing/immunity_composition.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/bell/immunity_composition.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/U/immunity_composition.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/perfect_matching_all/immunity_composition.png" alt="Drawing" style="width: 250px;"/> </td>
</tr></table>

#### Pathogen composition
<table><tr>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/decaying/pathogen_composition.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/increasing/pathogen_composition.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/bell/pathogen_composition.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/U/pathogen_composition.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/perfect_matching_all/pathogen_composition.png" alt="Drawing" style="width: 250px;"/> </td>
</tr></table>

#### Composition
<table><tr>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/decaying/weights_dist_composition.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/increasing/weights_dist_composition.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/bell/weights_dist_composition.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/U/weights_dist_composition.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/perfect_matching_all/weights_dist_composition.png" alt="Drawing" style="width: 250px;"/> </td>
</tr></table>

#### Compartments
<table><tr>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/decaying/reassortment_compartments.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/increasing/reassortment_compartments.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/bell/reassortment_compartments.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/U/reassortment_compartments.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/perfect_matching_all/reassortment_compartments.png" alt="Drawing" style="width: 250px;"/> </td>
</tr></table>

#### Infections
<table><tr>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/decaying/populations_plot.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/increasing/populations_plot.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/bell/populations_plot.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/U/populations_plot.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/perfect_matching_all/populations_plot.png" alt="Drawing" style="width: 250px;"/> </td>
</tr></table>

#### Cluster map
<table><tr>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/decaying/weights_dist_clustermap.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/increasing/weights_dist_clustermap.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/bell/weights_dist_clustermap.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/U/weights_dist_clustermap.png" alt="Drawing" style="width: 250px;"/> </td>
<td> <img src="immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on/perfect_matching_all/weights_dist_clustermap.png" alt="Drawing" style="width: 250px;"/> </td>
</tr></table>