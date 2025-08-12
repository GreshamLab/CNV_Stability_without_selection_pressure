

###Load libraries

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats
from matplotlib.gridspec import GridSpec

#--------------------------------------------

### Figure 2 

fig = plt.figure(figsize=(6,7),layout="constrained")

gs = GridSpec(4, 2, figure=fig)
ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[1,0])
ax3 = fig.add_subplot(gs[0:2,1])
ax4 = fig.add_subplot(gs[2:4,0:2])

#-------------

ax = ax1

df = pd.read_csv('Pairwise_Comps_III_Stats_Python.CSV')

Strain_names = ['2543_1','2543_2','2543_3']

color_array = ['navy','deepskyblue','blue']
linestyle_array = ['o:','s:','^:']

for j in range(len(Strain_names)):
    gen_names = ['Gen_0','Gen_2_30C','Gen_5_30C','Gen_8_30C','Gen_12_30C','Gen_15_30C'] 
    gen = np.array([0,2,5,8,12,15]) 
    val = np.NaN * gen
    for i in range(len(val)):
        val[i] = df[df['Sample']==Strain_names[j]][gen_names[i]]
    ax.plot(gen,(100-val), linestyle_array[j], color=color_array[j],ms=4)

ax.set_ylabel('Abundance of \n strain (percent)')
ax.set_ylim(15,55)


#-------------

ax = ax2

color_array = ['navy','deepskyblue','blue','navy','deepskyblue','blue','navy','deepskyblue','blue']
point_styles = ['o','s','^','o','s','^','o','s','^']

Strain ='DGY2543'

Strain_names = ['2543_1','2543_2','2543_3'] 
Gen_values = np.array([0,2,5,8,12,15,0,2,5,8,12,15,0,2,5,8,12,15])

val_combine_30 = np.array([])    
for k in range(len(Strain_names)):
    gen_names = ['Gen_0','Gen_2_30C','Gen_5_30C','Gen_8_30C','Gen_12_30C','Gen_15_30C']
    gen = np.array([0,2,5,8,12,15]) 
    val = np.NaN * gen
    Ln_ratio = np.NaN * gen
    Ln_adj = np.NaN * gen
    for i in range(len(val)):
        val[i] = df[df['Sample']==Strain_names[k]][gen_names[i]]
        Ln_ratio[i] = np.log((100-val[i])/val[i])
        Ln_adj[i] = Ln_ratio[i]-Ln_ratio[0]

    val_combine_30 = np.concatenate((val_combine_30,Ln_adj))

    ax.plot(gen,Ln_adj,point_styles[k],c=color_array[k],alpha=1.0,ms=4)

slope, intercept, r_value, p_value, std_err = stats.linregress(Gen_values, val_combine_30)

ax.plot(Gen_values, intercept + slope*Gen_values, c='royalblue')
ax.text(4,-0.3, "Slope = -0.065"'\n'"95% C.I. = -0.08, -0.05")

ax.set_xlabel('Generations')
ax.set_ylabel('ln (Strain/Ancestor)')

#-------------

ax = ax3

df1 = pd.read_csv('Fitness_vs_CNVsize.CSV')

ax.errorbar(df1['Extra_bp']/1000, df1['Slope']+1, yerr=2*df1['Std_err'], fmt='o',ms=4,
            color='black', ecolor='steelblue', capsize=3, lw=1.5)

ax.set_xlabel("CNV size (kilobases)")
ax.set_ylabel("Relative fitness")
ax.set_ylim(0.92, 1.11)
ax.set_xlim(10, 700)
ax.set_xscale('log')
ax.set_xticks([10,20,50,100,200,500],[10,20,50,100,200,500])

# Monte Carlo simulation to account for error bars
n_simulations = 10000
simulated_corrs = []
simulated_slope = []

x_arr = np.linspace(0,600,10)

for _ in range(n_simulations):
    simulated_slopes = np.random.normal(df1['Slope'], df1['Std_err'])  # Sample from normal dist using errors
    sim_corr, _ = stats.pearsonr(np.log10(df1['Extra_bp']), simulated_slopes)
    simulated_corrs.append(sim_corr)

# Compute mean and confidence intervals
corr_mean = np.mean(simulated_corrs)
lower, upper = np.percentile(simulated_corrs, [2.5, 97.5])

print(f"Mean Pearson correlation: {corr_mean:.3f}")
print(f"95% CI: ({lower:.3f}, {upper:.3f})")

ax.plot([0,1000],[1,1],':k')

ax.set_yticks([0.95,1.00,1.05,1.1])

#-------------

ax = ax4

df = pd.read_csv('Pairwise_Regression_Norm_All_30C_final.CSV')

mask1 = (df["Type"] == 'Rev')
mask2 = (df["Type"] == 'Rev_partial')

df1 = df[~mask1 & ~mask2]

ax.bar(x = df1['Strain'], height = df1['Slope']+1.,
      facecolor='gainsboro',edgecolor='black',width=0.5)


mask = (df1['P_value']<=0.05)
ax.plot(np.where(mask)[0], df1[mask]['Slope'] + 1. + 0.05,'*k')

# Add error bars
x_pos = np.arange(len(df1))  # Get the positions for the bars
ax.errorbar(x=x_pos, y=df1['Slope']+1., yerr=2*df1['Std_err'], fmt='none', ecolor='steelblue', capsize=3, lw=2)
ax.set_xticks(x_pos, df1['Strain'], rotation='vertical')

ax.set_ylim(0.850001,1.09999)
ax.set_xlabel('Strains', labelpad=28, fontsize=10.5)
ax.set_ylabel('Relative fitness')

ax.plot([-0.5,np.max(x_pos)+0.5],[1.,1],'--',color='brown')
ax.set_xlim(-0.55,np.max(x_pos)+0.55)

# Segmental CNVs 
x_start_seg = 0
x_end_seg = 10
y_bracket = 0.775  
ax.plot([x_start_seg, x_start_seg, x_end_seg, x_end_seg], 
        [y_bracket+0.02, y_bracket, y_bracket, y_bracket+0.02], 
        lw=1.5, color='grey', clip_on=False)
ax.text((x_start_seg + x_end_seg)/2, y_bracket + 0.015, 'Segmental CNVs', 
        ha='center', va='top', fontsize=8.5)

# Aneuploid 
x_start_aneu = 11
x_end_aneu = 14
ax.plot([x_start_aneu, x_start_aneu, x_end_aneu, x_end_aneu], 
        [y_bracket+0.01, y_bracket, y_bracket, y_bracket+0.01], 
        lw=1.5, color='grey', clip_on=False)
ax.text((x_start_aneu + x_end_aneu)/2, y_bracket + 0.015, 'Aneuploidies', 
        ha='center', va='top', fontsize=8.5)


#plt.savefig('Fig2.pdf',bbox_inches='tight',pad_inches=0.5)
#plt.savefig('Fig2.png',pad_inches=0.05)

plt.show()

#--------------------------------------------


