import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd

path_save = 'tentative_figs/fig3_immunity_effect_over_evolution/immunity_effect_over_evolution_weights_dist_simple/not_grouped_recomb_on'

df = pd.DataFrame()
df['pos'] = [0,1,2]
df['weight'] = [1,1,1]
fig, ax = plt.subplots(1, 1, figsize=(5,5))
sns.barplot(ax=ax, data=df,x='pos', y='weight', label='perfect_matching_all')
#eliminate x axis
ax.set_xlabel('')
ax.set_ylim(0,1.01)
# ax.set_xticks([])
ax.set_ylabel('Weight')
fig.savefig(os.path.join(path_save,'weights_dist_perfect_matching.png'))