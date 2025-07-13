

###Load libraries

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#-----------------------------------

###Figure 6: Empirical fitness of revertants and CNV parent strains, alongwith fitness estimates from simulation-based inference

df = pd.read_csv('Pairwise_Regression_Norm_All_30C_final.CSV')

mask1 = (df["Type"] == 'CNV')
df1 = df[~mask1]

ax = df1.plot(x='Strain',              
             y='Slope_plus_one', 
             kind='bar',                    
             figsize=(12,6),          
             ylabel='Relative fitness',           
             rot=90,
             color='grey',
             edgecolor='k',
             label=None)  
            
c=(np.where(df1['Type']=='Rev')[0])
print(c)
bars = ax.patches
for i in range(len(c)):
    bars[c[i]].set_facecolor('gainsboro');
    
    
c=(np.where(df1['Type']=='Rev_partial')[0])
print(c)
bars = ax.patches
for i in range(len(c)):
    bars[c[i]].set_facecolor('gainsboro')
    bars[c[i]].set_hatch('////')

# Add error bars 
x_pos = np.arange(len(df1))  # Get the positions for the bars
ax.errorbar(x=x_pos, y=df1['Slope_plus_one'], yerr=2*df1['Std_err'], fmt='none', ecolor='steelblue', capsize=3, lw=2)
ax.legend('')

ax.set_xticks([0,11,18,23,32],['G_4','GM_12','GM_13','GM_14','GM_15'])


ax.annotate('', xy=(0.075, 0.11), xytext=(0.245, 0.11), xycoords='figure fraction',
                  arrowprops=dict(arrowstyle='-', color='k', linewidth=1.2, linestyle='-'))

ax.annotate('', xy=(0.265, 0.11), xytext=(0.37, 0.11), xycoords='figure fraction',
                  arrowprops=dict(arrowstyle='-', color='k', linewidth=1.2, linestyle='-'))

ax.annotate('', xy=(0.39, 0.11), xytext=(0.455, 0.11), xycoords='figure fraction',
                  arrowprops=dict(arrowstyle='-', color='k', linewidth=1.2, linestyle='-'))

ax.annotate('', xy=(0.48, 0.11), xytext=(0.615, 0.11), xycoords='figure fraction',
                  arrowprops=dict(arrowstyle='-', color='k', linewidth=1.2, linestyle='-'))

ax.annotate('', xy=(0.64, 0.11), xytext=(0.825, 0.11), xycoords='figure fraction',
                  arrowprops=dict(arrowstyle='-', color='k', linewidth=1.2, linestyle='-'))

ax.text(4,0.85,'G_4_Rev')

ax.text(12,0.85,'GM_12_Rev')

ax.text(18.5,0.85,'GM_13_Rev')

ax.text(25.5,0.85,'GM_14_Rev')

ax.text(36,0.85,'GM_15_Rev')

ax.plot([10.5,10.5],[0.87,1.13],':k')
ax.plot([17.5,17.5],[0.87,1.13],':k')
ax.plot([22.5,22.5],[0.87,1.13],':k')
ax.plot([31.5,31.5],[0.87,1.13],':k')
ax.plot([-0.5,np.max(x_pos)+0.5],[1.,1],'--',color='brown')

ax.set_ylim(0.87,1.13)

# Maximum and minimum estimates for each CNV strain, obtained from simulation-based inference
highlight_estimates = {
    'G_4': (0.9798, 0.9694),
    'GM_12': (0.9476757, 0.92209492),
    'GM_13': (0.9534311, 0.92233656),
    'GM_14': (0.92137516, 0.90289836),
    'GM_15': (0.9444589, 0.91573394)
}

# Reset index to get positional index corresponding to plotting
df1 = df1.reset_index(drop=True)

for strain, (ymax, ymin) in highlight_estimates.items():
    if strain in df1['Strain'].values:
        x = df1[df1['Strain'] == strain].index[0]
        ax.fill_between([x - 0.3, x + 0.3], ymin, ymax, color='cyan', alpha=0.6)


#plt.savefig('Fig6.png',bbox_inches='tight',pad_inches = 0.05,dpi=300)
#plt.savefig('Fig6.pdf',bbox_inches='tight',pad_inches = 0.05,dpi=300)

plt.show()

#-----------------------------------

