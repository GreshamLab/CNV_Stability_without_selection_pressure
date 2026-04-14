

###Load libraries

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable

#--------------------------------------------

### Figure 2 

fig = plt.figure(figsize=(6.5,7),layout="constrained")

gs = GridSpec(4, 2, figure=fig)
ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[1,0])
ax3 = fig.add_subplot(gs[0:2,1])
ax4 = fig.add_subplot(gs[2:4,0:2])

#-------------

ax = ax1

df = pd.read_csv('Pairwise_III_edited_new.CSV')

Strain_names = ['2750_1','2750_2','2750_3']

color_array = ['navy','deepskyblue','blue']
linestyle_array = ['^:','o:','s:']

for j in range(len(Strain_names)):
    gen_names = ['Gen_0','Gen_2_30C','Gen_5_30C','Gen_8_30C','Gen_12_30C','Gen_15_30C'] 
    gen = np.array([0,2,5,8,12,15]) 
    val = np.NaN * gen
    for i in range(len(val)):
        val[i] = df[df['Sample']==Strain_names[j]][gen_names[i]]
    ax.plot(gen,(100-val), linestyle_array[j], color=color_array[j],ms=4)

ax.set_ylabel('Abundance of \n strain (percent)')
ax.set_ylim(15,80)
ax.text(4.5,70,'3 replicates of GM_13')


#-------------

ax = ax2

color_array = ['navy','deepskyblue','blue','navy','deepskyblue','blue','navy','deepskyblue','blue']
point_styles = ['^','o','s','^','o','s','^','o','s']

Strain ='DGY2750'

Strain_names = ['2750_1','2750_2','2750_3'] 
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
print(slope)
print(std_err)

ax.plot(Gen_values, intercept + slope*Gen_values, c='royalblue')
ax.text(5.3,-0.4, "Slope = -0.086"'\n'"95% C.I.: (-0.05, -0.12)")

ax.set_xlabel('Generations')
ax.set_ylabel('ln (Strain/Ancestor)')

#-------------

ax = ax3

df1 = pd.read_csv('GaschFitnessCost_vs_Size.CSV')

ax.errorbar(df1['CNV_cost'], df1['Slope']+1, yerr=2*df1['Std_err'], fmt='o',ms=4,
            color='black', ecolor='grey', capsize=3, lw=1.5)


ax.set_xlabel("             Calculated CNV cost")
ax.set_ylabel("Relative fitness")
ax.set_ylim(0.9, 1.11)
ax.set_xlim(-82,1)
ax.set_yticks([0.95,1.00,1.05,1.1])


#Monte Carlo simulation to account for error bars
n_simulations = 10000
simulated_corrs = []
simulated_slope = []

x_arr = np.linspace(0,600,10)

for _ in range(n_simulations):
    simulated_slopes = np.random.normal(df1['Slope'], df1['Std_err']) 
    sim_corr, _ = stats.pearsonr(df1['CNV_cost'], simulated_slopes)
    simulated_corrs.append(sim_corr)

#Monte Carlo regression: compute slope + intercept for each simulation 
sim_slopes_mc = []
sim_intercepts_mc = []

for _ in range(n_simulations):
    y_sim = np.random.normal(df1['Slope']+1, df1['Std_err'])
    slope_mc, intercept_mc = np.polyfit(df1['CNV_cost'], y_sim, 1)
    sim_slopes_mc.append(slope_mc)
    sim_intercepts_mc.append(intercept_mc)

#Median regression line based on MC fits
median_slope = np.median(sim_slopes_mc)
median_intercept = np.median(sim_intercepts_mc)

x_line = np.linspace(df1['CNV_cost'].min(), df1['CNV_cost'].max(), 200)
y_line = median_slope * x_line + median_intercept

ax.plot(x_line, y_line, '--', color='black', lw=1.5, zorder=1)

#Linear regression confidence interval 

y_all = np.zeros((n_simulations, len(x_line)))
y_high = np.zeros_like(x_line)
y_low = np.zeros_like(x_line)
for i in range(n_simulations):
    y_all[i,:] = sim_slopes_mc[i] * x_line + sim_intercepts_mc[i]

for j in range(len(x_line)):
    y_high[j] = np.percentile(y_all[:,j], 97.5)
    y_low[j] = np.percentile(y_all[:,j], 2.5)
ax.fill_between(x_line, y_high, y_low, color='grey', alpha=0.3, zorder=-1)

#Compute mean and confidence intervals
corr_mean = np.mean(simulated_corrs)
lower, upper = np.percentile(simulated_corrs, [2.5, 97.5])

#print(f"Mean Pearson correlation: {corr_mean:.3f}")
#print(f"95% CI: ({lower:.3f}, {upper:.3f})")

ax.text(-80, 1.085, 'Mean Pearson correlation: 0.385\n95% C.I.: (0.153, 0.609)')

ax.annotate('', xy=(0.10, -0.135), xycoords='axes fraction', xytext=(0.3, -0.135), 
            arrowprops=dict(arrowstyle="->", color='k'))

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

#plt.savefig('Fig2_revision2.pdf',bbox_inches='tight',pad_inches=0.5)
#plt.savefig('Fig2.png',pad_inches=0.05, dpi=600)

plt.show()

#--------------------------------------------


