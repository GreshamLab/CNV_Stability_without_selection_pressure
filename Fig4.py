

###Load libraries

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import statistics
from scipy import interpolate
from scipy.interpolate import CubicSpline
from itertools import combinations
from scipy.stats import ttest_ind
from matplotlib.gridspec import GridSpec
from scipy.optimize import root_scalar

#-----------------------------------

###Read in flow cytometry data from LTEE

df_2 = pd.read_csv('2023_12_Percentage_gating_LTEE2_edit.csv')

#-----------------------------------

### Figure 4

fig = plt.figure(figsize=(6,6),layout="constrained")

gs = GridSpec(2, 2, figure=fig)
ax0 = fig.add_subplot(gs[0,0])
ax1 = fig.add_subplot(gs[0,1])
ax2 = fig.add_subplot(gs[1,0])
ax3 = fig.add_subplot(gs[1,1])

ax = ax0

df2 = df_2.loc[df_2['Strain'] == 'DGY1886_1']
data = df2[['Generation', 'Multicopy_mCitrine']].to_numpy()
    
x= data[:,0]
y= data[:,1]

# Create the cubic spline interpolation
cs = CubicSpline(x,y)
# Generate x values for the smooth curve
x_smooth = np.linspace(x.min(), x.max(), 100)
# Evaluate the spline at the smooth x values
y_smooth = cs(x_smooth)

x_75 = 220
for x_val_1 in x_smooth:
    if cs(x_val_1) <= 75 and cs(x_val_1) > cs(x_75):
        x_75 = x_val_1

x_50 = 220
for x_val_2 in x_smooth:
    if cs(x_val_2) <= 50 and cs(x_val_2) > cs(x_50):
        x_50 = x_val_2

x_25 = 220
for x_val_3 in x_smooth:
    if cs(x_val_3) <= 25 and cs(x_val_3) > cs(x_25):
        x_25 = x_val_3

ax.plot(x, y, '.', color='mediumblue') 
ax.plot(x_smooth, y_smooth, color='royalblue')

# Plot dotted lines at the intersection point
ax.plot([x_75, x_75], [-0, 75], ':',color='k', linewidth=1)  # Vertical line to x-axis
ax.plot([x_smooth[0], x_75], [75, 75], ':',color='k', linewidth=1)  # Horizontal line to curve

# Function to solve for x when cs(x) = 50
def find_x_at_y(target_y, cs, x_range):
    func = lambda x: cs(x) - target_y
    result = root_scalar(func, bracket=[x_range.min(), x_range.max()])
    return result.root if result.converged else None
x_50 = find_x_at_y(50, cs, x_smooth)

if x_50 is not None:
    ax.plot([x_smooth[0], x_50], [50, 50], '--',color='k', linewidth=1)  
    ax.plot([x_50, x_50], [0, 50], '--',color='k', linewidth=1)  

# Function to solve for x when cs(x) = 25
def find_x_at_y(target_y, cs, x_range):
    func = lambda x: cs(x) - target_y
    result = root_scalar(func, bracket=[x_range.min(), x_range.max()])
    return result.root if result.converged else None
x_25 = find_x_at_y(25, cs, x_smooth)

if x_25 is not None:
    ax.plot([x_smooth[0], x_25], [25, 25], '--',color='black', linewidth=1.5) 
    ax.plot([x_25, x_25], [0, 25], '--',color='black', linewidth=1.5)

ax.text(10, 77, 'Early phase', color='royalblue')
ax.text(10, 52, 'Middle phase', color='navy')
ax.text(10, 27, 'Late phase', color='k')

ax.set_xlabel('Generations')
ax.set_ylabel('Percentage of cells with CNV')
ax.set_ylim(0, 103)
ax.set_xlim(0, 225)
ax.set_yticks([0, 25, 50, 75, 100])
ax.set_xticks([0, 75, 150, 220])

#-------------------

##Data
#Early
G_5 = [106, 112, 116, 151]
GM_Chr11 = [38, 44, 53, 66, 68, 86, 99, 103, 134, 167, 202, 300, 300] 
GM_Chr14 = [7, 10, 10, 12, 16, 16, 18, 18, 18, 23, 58, 62, 64]

#Middle
G_5_m = [158, 149, 206, 204]
GM_Chr11_m = [58, 77, 88, 143, 151, 169, 193, 206, 300, 300, 300, 300, 300]
GM_Chr14_m = [18, 31, 36, 36, 36, 38, 38, 40, 40, 40, 68, 75, 75]

#Late
G_5_l = [212, 208, 300, 300]
GM_Chr11_l = [140, 149, 162, 188, 300, 300, 300, 300, 300, 300, 300, 300, 300]
GM_Chr14_l = [36, 44, 49, 49, 49, 49, 51, 51,  51, 53, 79, 86, 86]

## Pairwise t-tests
groups = [G_5, GM_Chr11, GM_Chr14]
groups_m = [G_5_m, GM_Chr11_m, GM_Chr14_m]
groups_l = [G_5_l, GM_Chr11_l, GM_Chr14_l]

labels = ['Chr11_partial', 'Chr11_aneuploid', 'Chr14_aneuploid']
labels_m = ['Chr11_partial_m', 'Chr11_aneuploid_m', 'Chr14_aneuploid_m']
labels_l = ['Chr11_partial_l', 'Chr11_aneuploid_l', 'Chr14_aneuploid_l']

#Early
p_values = {}
significant_pairs = []
for (i, j) in combinations(range(len(groups)), 2):
    stat, p = ttest_ind(groups[i], groups[j], equal_var=False)  # Welch’s t-test
    p_values[(labels[i], labels[j])] = p
    if p < 0.05:
        significant_pairs.append((i, j))
    
#Middle
p_values = {}
significant_pairs = []
for (i, j) in combinations(range(len(groups_m)), 2):
    stat, p = ttest_ind(groups_m[i], groups_m[j], equal_var=False)  
    p_values[(labels_m[i], labels_m[j])] = p
    if p < 0.05:
        significant_pairs.append((i, j))
    
#Late
p_values = {}
significant_pairs = []
for (i, j) in combinations(range(len(groups_l)), 2):
    stat, p = ttest_ind(groups_l[i], groups_l[j], equal_var=False)  
    p_values[(labels_l[i], labels_l[j])] = p
    if p < 0.05:
        significant_pairs.append((i, j))

## Box plot
#Early
filtered_data = [[val for val in group if val != 300] for group in [G_5, GM_Chr11, GM_Chr14]]
excluded_points = [(i + 1, val) for i, group in enumerate([G_5, GM_Chr11, GM_Chr14]) for val in group if val == 300]

box_colors = ['springgreen', 'greenyellow', 'orange']  
boxplots = ax1.boxplot(filtered_data, labels=['Chr XI', 'Chr XI', 'Chr XIV'], patch_artist=True, showfliers=False)

i = 1
for y_val in filtered_data:
    ax1.plot(np.array(y_val)*0.+i,y_val,'.',mfc='none',mec='k')
    i = i+1

for patch, color in zip(boxplots['boxes'], box_colors):
    patch.set_facecolor(color)
for median in boxplots['medians']:
        median.set_color('black')
        median.set_linewidth(1.5)

#Writing significant p-values below boxplots    
x2, x3 = 2, 3  #Positions of 'Chr11_aneuploid' and 'Chr14_aneuploid'

y, h = -30, -10 
ax1.plot([x2, x2, x3, x3], [y, y + h, y + h, y], color='brown', linewidth=1.2)
ax1.text(0.5*(x2+x3), y-5, 'p = 0.0013', ha='center', va='bottom')

x1, x3 = 1, 3  #Positions of 'Chr11_partial' and 'Chr14_aneuploid'
y, h = -100, -10
ax1.plot([x1, x1, x3, x3], [y, y + h, y + h, y], color='brown', linewidth=1.2)
ax1.text(0.5*(x1+x3), y-5, 'p = 0.0004', ha='center', va='bottom')

for x_pos, y_val in excluded_points:
    ax1.plot(x_pos + np.random.uniform(-0.15,0.15), y_val, marker='^', color='g', markersize=6, alpha=0.4) 

ax1.set_title('Early phase', color='royalblue')


#Middle
filtered_data_m = [[val for val in group if val != 300] for group in [G_5_m, GM_Chr11_m, GM_Chr14_m]]
excluded_points_m = [(i + 1, val) for i, group in enumerate([G_5_m, GM_Chr11_m, GM_Chr14_m]) for val in group if val == 300]

boxplots = ax2.boxplot(filtered_data_m, labels=['Chr XI', 'Chr XI', 'Chr XIV'], patch_artist=True, showfliers=False)

i = 1
for y_val in filtered_data_m:
    ax2.plot(np.array(y_val)*0.+i,y_val,'.',mfc='none',mec='k')
    i = i+1
    
for patch, color in zip(boxplots['boxes'], box_colors):
    patch.set_facecolor(color)
for median in boxplots['medians']:
        median.set_color('black')
        median.set_linewidth(1.5)
        
x2, x3 = 2, 3 
y, h = -30, -10 
ax2.plot([x2, x2, x3, x3], [y, y + h, y + h, y], color='brown', linewidth=1.2)
ax2.text(0.5*(x2+x3), y-5, 'p = 0.0001', ha='center', va='bottom')

x1, x3 = 1, 3 
y, h = -100, -10
ax2.plot([x1, x1, x3, x3], [y, y + h, y + h, y], color='brown', linewidth=1.2)
ax2.text(0.5*(x1+x3), y-5, 'p = 0.0016', ha='center', va='bottom')

for x_pos, y_val in excluded_points_m:
    ax2.plot(x_pos + np.random.uniform(-0.15,0.15), y_val, marker='^', color='g', markersize=6, alpha=0.4) 

ax2.set_title('Middle phase', color='navy')

#Late
filtered_data_l = [[val for val in group if val != 300] for group in [G_5_l, GM_Chr11_l, GM_Chr14_l]]
excluded_points_l = [(i + 1, val) for i, group in enumerate([G_5_l, GM_Chr11_l, GM_Chr14_l]) for val in group if val == 300]

boxplots = ax3.boxplot(filtered_data_l, labels=['Chr XI', 'Chr XI', 'Chr XIV'], patch_artist=True, showfliers=False)

i = 1
for y_val in filtered_data_l:
    ax3.plot(np.array(y_val)*0.+i,y_val,'.',mfc='none',mec='k')
    i = i+1
    
for patch, color in zip(boxplots['boxes'], box_colors):
    patch.set_facecolor(color)
for median in boxplots['medians']:
        median.set_color('black')
        median.set_linewidth(1.5)
        
x2, x3 = 2, 3  
y, h = -30, -10 
ax3.plot([x2, x2, x3, x3], [y, y + h, y + h, y], color='brown', linewidth=1.2)
ax3.text(0.5*(x2+x3), y-5, 'p = 0.0000', ha='center', va='bottom')

x1, x3 = 1, 3  
y, h = -100, -10
ax3.plot([x1, x1, x3, x3], [y, y + h, y + h, y], color='brown', linewidth=1.2)
ax3.text(0.5*(x1+x3), y-5, 'p = 0.0041', ha='center', va='bottom')

for x_pos, y_val in excluded_points_l:
    ax3.plot(x_pos + np.random.uniform(-0.15,0.15), y_val, marker='^', color='g', markersize=6, alpha=0.5) 

ax3.set_title('Late phase')

#Common for all 3 plots
for ax in [ax1,ax2,ax3]:
    ax.set_ylim(-150,320)
    ax.set_xlim(0.8,3.2)
    ax.plot([0,4],[0,0],':k')
    ax.plot([0,4],[220,220],':k')
    ax.set_yticks(ticks=[-100,0,50,100,150,200,250,300],labels=['', '0','50','100','150','200','250','300'])
    plt.setp(ax.get_xticklabels(), rotation=0)
    ax.yaxis.set_tick_params(labelleft=True)
    ax.set_yticks(ticks=[0,75,150,220],labels=['0','75','150','220'])
    ax.set_ylabel('Generations')

    # Segmental and aneuploidy labels 
    x_start_seg = 0.8
    x_end_seg = 1.2
    y_bracket = -220  
    ax.plot([x_start_seg, x_start_seg, x_end_seg, x_end_seg], 
        [y_bracket+10, y_bracket, y_bracket, y_bracket+10], 
        lw=1.5, color='k', clip_on=False)
    ax.text(1, y_bracket-8, 'Segmental', 
        ha='center', va='top', fontsize=10)
    
    x_start_aneu = 2
    x_end_aneu = 3
    y_bracket = -220  
    ax.plot([x_start_aneu, x_start_aneu, x_end_aneu, x_end_aneu], 
        [y_bracket+10, y_bracket, y_bracket, y_bracket+10], 
        lw=1.5, color='k', clip_on=False)
    ax.text((x_start_aneu + x_end_aneu)/2, y_bracket-8, 'Aneuploid', 
        ha='center', va='top', fontsize=10)

plt.tight_layout()

# Define datasets for median calculation
datasets = {
    "G_5_early": G_5, "G_5_middle": G_5_m, "G_5_late": G_5_l,
    "GM_Chr11_early": GM_Chr11, "GM_Chr11_middle": GM_Chr11_m, "GM_Chr11_late": GM_Chr11_l,
    "GM_Chr14_early": GM_Chr14, "GM_Chr14_middle": GM_Chr14_m, "GM_Chr14_late": GM_Chr14_l
}

# Compute and print medians excluding 300
medians = {name: statistics.median([val for val in data if val != 300]) for name, data in datasets.items()}
print(medians)

#plt.savefig('Fig4.png', pad_inches= 0.5, dpi=300)
#plt.savefig('Fig4.pdf', pad_inches= 0.5, dpi=300)

plt.show()

#-----------------------------------
