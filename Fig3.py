

###Load required packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import statistics
import matplotlib.patches as patches

#-----------------------------------

###Read in flow cytometry data from LTEE

df_1 = pd.read_csv('2024_11_Percentage_gating_LTEE3_YPD_20241205_edited.csv')
df_2 = pd.read_csv('2023_12_Percentage_gating_LTEE2_edit.csv')
df_0 = pd.read_csv('2018_Percentage_CNV_acquire_Lauer.csv')
df_02 = pd.read_csv('2018_Percentage_CNV_acquire_Rahman.csv')
df_03 = pd.read_csv('2018_Percentage_CNV_acquire_Avecilla.csv')

#-----------------------------------

###First panel

fig, axs = plt.subplots(1,4, sharex=True, sharey=True, figsize=(8,2.5))
fig.subplots_adjust(wspace=0.1,hspace=0.3)

gb_cols = ['limegreen', 'green', 'darkturquoise']
ctrl_cols = ['darkgrey', 'black','grey']

#-----------

jj=0

strain_names = ['DGY1741_1','DGY1741_2','DGY1741_3']

shift_values = np.linspace(-2, 2, len(strain_names))  # Generate small shifts for each replicate

for j in range(len(strain_names)):
    df2 = df_1.loc[df_1['Strain'] == strain_names[j]]
    data = df2[['Generation','Multicopy_mCitrine']].to_numpy()
    axs[jj].plot( data[:,0],data[:,1] + shift_values[j],'.-',color=gb_cols[j],lw=2)

    
axs[jj].text(40,117,'G_1')
axs[jj].set_ylim(-5,110)

##Inset
inset_ax = axs[jj].inset_axes([0.18, 0.14, 0.4, 0.35])  # [x, y, width, height] relative to axs[0]
strain_names = ['DGY1741']

for j in range(len(strain_names)):
    df2 = df_0.loc[df_0['Strain'] == strain_names[j]]
    data = df2[['Generation', 'Multicopy_mCitrine']].to_numpy()
    inset_ax.plot(data[:, 0], data[:, 1], '.-', lw=2, markersize=0.5)
    inset_ax.set_facecolor('khaki')
    
#Mark the generation at which the strain was isolated
    gen_iso = df2.loc[df2['Generation'] == 244]
    inset_ax.scatter(gen_iso['Generation'], gen_iso['Multicopy_mCitrine'], marker='o', facecolor = 'k', edgecolor='k', s=10, zorder=10)


inset_ax.set_ylim(-5, 110)
inset_ax.set_yticks([0,100], [0,100], fontsize=6.5)
inset_ax.set_xticks([0,250], [0,250], fontsize=6.5)

#-----------

jj=1

strain_names = ['DGY1743_1','DGY1743_2','DGY1743_3']

shift_values = np.linspace(-2, 2, len(strain_names)) 

for j in range(len(strain_names)):
    df2 = df_1.loc[df_1['Strain'] == strain_names[j]]
    data = df2[['Generation','Multicopy_mCitrine']].to_numpy()
    axs[jj].plot( data[:,0],data[:,1] + shift_values[j],'.-',color=gb_cols[j],lw=2)

axs[jj].text(40,117,'G_2')
axs[jj].set_ylim(-5,110)

##Inset
inset_ax = axs[jj].inset_axes([0.18, 0.14, 0.4, 0.35]) 
strain_names = ['DGY1743']

for j in range(len(strain_names)):
    df2 = df_0.loc[df_0['Strain'] == strain_names[j]]
    data = df2[['Generation', 'Multicopy_mCitrine']].to_numpy()
    inset_ax.plot(data[:, 0], data[:, 1], '.-', lw=2, markersize=0.5)
    inset_ax.set_facecolor('khaki')
    
#Mark the generation at which the strain was isolated
    gen_iso = df2.loc[df2['Generation'] == 244]
    inset_ax.scatter(gen_iso['Generation'], gen_iso['Multicopy_mCitrine'], marker='o', facecolor = 'k', edgecolor='k', s=10, zorder=10)


inset_ax.set_ylim(-5, 110)
inset_ax.set_yticks([0,100], [0,100], fontsize=6.5)
inset_ax.set_xticks([0,250], [0,250], fontsize=6.5)

#-----------

jj=2

strain_names = ['DGY1744_1','DGY1744_2','DGY1744_3']

shift_values = np.linspace(-2, 2, len(strain_names)) 

for j in range(len(strain_names)):
    df2 = df_1.loc[df_1['Strain'] == strain_names[j]]
    data = df2[['Generation','Multicopy_mCitrine']].to_numpy()
    axs[jj].plot( data[:,0],data[:,1] + shift_values[j],'.-',color=gb_cols[j],lw=2)

    
axs[jj].text(40,117,'G_3')
axs[jj].set_ylim(-5,110)

##Inset
inset_ax = axs[jj].inset_axes([0.18, 0.14, 0.4, 0.35])  
strain_names = ['DGY1744']

for j in range(len(strain_names)):
    df2 = df_0.loc[df_0['Strain'] == strain_names[j]]
    data = df2[['Generation', 'Multicopy_mCitrine']].to_numpy()
    inset_ax.plot(data[:, 0], data[:, 1], '.-', lw=2, markersize=0.5)
    inset_ax.set_facecolor('khaki')
    
#Mark the generation at which the strain was isolated
    gen_iso = df2.loc[df2['Generation'] == 244]
    inset_ax.scatter(gen_iso['Generation'], gen_iso['Multicopy_mCitrine'], marker='o', facecolor = 'k', edgecolor='k', s=10, zorder=10)


inset_ax.set_ylim(-5, 110)
inset_ax.set_yticks([0,100], [0,100], fontsize=6.5)
inset_ax.set_xticks([0,250], [0,250], fontsize=6.5)

#-----------

jj=3


strain_names = ['DGY1315_1','DGY1315_2','DGY1315_3']

shift_values = np.linspace(-2, 2, len(strain_names)) 

for j in range(len(strain_names)):
    df2 = df_1.loc[df_1['Strain'] == strain_names[j]]
    data = df2[['Generation','Multicopy_mCitrine']].to_numpy()
    axs[jj].plot( data[:,0],data[:,1] + shift_values[j],'.-',color=ctrl_cols[j],lw=2)
    
strain_names = ['DGY500_1','DGY500_2','DGY500_3']

shift_values = np.linspace(-2, 2, len(strain_names))  

for j in range(len(strain_names)):
    df2 = df_1.loc[df_1['Strain'] == strain_names[j]]
    data = df2[['Generation','Multicopy_mCitrine']].to_numpy()
    axs[jj].plot( data[:,0],data[:,1] + shift_values[j],'.-',color=ctrl_cols[j],lw=2)

axs[jj].text(4,85,'Multi-copy control', fontsize=9)
axs[jj].text(4,15,'One-copy control', fontsize=9)
axs[jj].text(17,117,'Control strains')
axs[jj].set_ylim(-5,110)

#-----------

axs[0].set_ylabel('Percentage of cells with CNV')
axs[1].set_xlabel('                                       Generations')

plt.tight_layout()

#plt.savefig('Fig_3A.png', pad_inches= 0.05, dpi=300)
#plt.savefig('Fig_3A.pdf', pad_inches= 0.05, dpi=300)

plt.show()

#-----------------------------------

###Second panel

fig, axs = plt.subplots(2,4, sharex=True, sharey=True, figsize=(8,4))
fig.subplots_adjust(wspace=0.1,hspace=0.3)

gb_cols = ['royalblue', 'limegreen', 'green', 'darkturquoise' ]
r_cols = ['peru', 'brown', 'orange', 'orangered']
p_cols = ['deeppink', 'mediumvioletred', 'pink', 'violet']

#-----------

ii=0
jj=0

strain_names = ['DGY1886_1','DGY1886_2','DGY1886_3','DGY1886_4']

for j in range(len(strain_names)):
    df2 = df_2.loc[df_2['Strain'] == strain_names[j]]
    data = df2[['Generation','Multicopy_mCitrine']].to_numpy()
    axs[ii,jj].plot( data[:,0],data[:,1] ,'.-',color=gb_cols[j],lw=2)
    
axs[ii,jj].text(80,115,'G_4')
axs[ii,jj].set_ylim(-5,110)

##Inset
inset_ax = axs[ii,jj].inset_axes([0.18, 0.14, 0.4, 0.35]) 
strain_names = ['DGY1883']

for j in range(len(strain_names)):
    df2 = df_03.loc[(df_03['Strain'] == strain_names[j]) & (df_03['Generation'] <= 150)]
    data = df2[['Generation', 'Multicopy_mCitrine']].to_numpy()
    inset_ax.plot(data[:, 0], data[:, 1], '.-', lw=2, markersize=0.5)
    inset_ax.set_facecolor('khaki')
    
#Mark the generation at which the strain was isolated
    gen_iso = df2.loc[df2['Generation'] == 70]
    inset_ax.scatter(gen_iso['Generation'], gen_iso['Multicopy_mCitrine'], marker='o', facecolor = 'k', edgecolor='k', s=10, zorder=10)
    
inset_ax.set_ylim(-5, 110)
inset_ax.set_xlim(-5, 160)
inset_ax.set_yticks([0,100], [0,100], fontsize=6.5)
inset_ax.set_xticks([0,150], [0,150], fontsize=6.5)

#-----------

ii=0
jj=1

strain_names = ['DGY1917_1','DGY1917_2','DGY1917_3','DGY1917_4']

shift_values = np.linspace(-3.5, 3.5, len(strain_names)) 

for j in range(len(strain_names)):
    df2 = df_2.loc[df_2['Strain'] == strain_names[j]]
    data = df2[['Generation','Multicopy_mCitrine']].to_numpy()
    axs[ii, jj].plot(data[:, 0], data[:, 1] + shift_values[j], '.-', color=gb_cols[j], lw=2)

axs[ii, jj].text(80, 115, 'G_5')
axs[ii, jj].set_ylim(-5, 110)


##Inset
inset_ax = axs[ii,jj].inset_axes([0.18, 0.14, 0.4, 0.35]) 
strain_names = ['DGY1883']

for j in range(len(strain_names)):
    df2 = df_03.loc[(df_03['Strain'] == strain_names[j]) & (df_03['Generation'] <= 150)]
    data = df2[['Generation', 'Multicopy_mCitrine']].to_numpy()
    inset_ax.plot(data[:, 0], data[:, 1], '.-', lw=2, markersize=0.5)
    inset_ax.set_facecolor('khaki')
    
#Mark the generation at which the strain was isolated
    gen_iso = df2.loc[df2['Generation'] == 70]
    inset_ax.scatter(gen_iso['Generation'], gen_iso['Multicopy_mCitrine'], marker='o', facecolor = 'k', edgecolor='k', s=10, zorder=10)
    
inset_ax.set_ylim(-5, 110)
inset_ax.set_xlim(-5, 160)
inset_ax.set_yticks([0,100], [0,100], fontsize=6.5)
inset_ax.set_xticks([0,150], [0,150], fontsize=6.5)

#-----------

ii=0
jj=2

strain_names = ['DGY1883_1','DGY1883_2','DGY1883_3','DGY1883_4']

shift_values = np.linspace(-3.5, 3.5, len(strain_names)) 

for j in range(len(strain_names)):
    df2 = df_2.loc[df_2['Strain'] == strain_names[j]]
    data = df2[['Generation','Multicopy_mCitrine']].to_numpy()
    axs[ii, jj].plot(data[:, 0], data[:, 1] + shift_values[j], '.-', color=gb_cols[j], lw=2)
    
axs[ii,jj].text(80,115,'G_6')
axs[ii,jj].set_ylim(-5,110)

##Inset
inset_ax = axs[ii,jj].inset_axes([0.18, 0.14, 0.4, 0.35]) 
strain_names = ['DGY1883']

for j in range(len(strain_names)):
    df2 = df_03.loc[(df_03['Strain'] == strain_names[j]) & (df_03['Generation'] <= 150)]
    data = df2[['Generation', 'Multicopy_mCitrine']].to_numpy()
    inset_ax.plot(data[:, 0], data[:, 1], '.-', lw=2, markersize=0.5)
    inset_ax.set_facecolor('khaki')
    
#Mark the generation at which the strain was isolated
    gen_iso = df2.loc[df2['Generation'] == 70]
    inset_ax.scatter(gen_iso['Generation'], gen_iso['Multicopy_mCitrine'], marker='o', facecolor = 'k', edgecolor='k', s=10, zorder=10)
    
inset_ax.set_ylim(-5, 110)
inset_ax.set_xlim(-5, 160)
inset_ax.set_yticks([0,100], [0,100], fontsize=6.5)
inset_ax.set_xticks([0,150], [0,150], fontsize=6.5)

#-----------

ii=0
jj=3

strain_names = ['DGY1988_1','DGY1988_2','DGY1988_3','DGY1988_4']

shift_values = np.linspace(-3.5, 3.5, len(strain_names))  

for j in range(len(strain_names)):
    df2 = df_2.loc[df_2['Strain'] == strain_names[j]]
    data = df2[['Generation','Multicopy_mCitrine']].to_numpy()
    axs[ii, jj].plot(data[:, 0], data[:, 1] + shift_values[j], '.-', color=gb_cols[j], lw=2)
    
axs[ii,jj].text(80,115,'G_7')
axs[ii,jj].set_ylim(-5,110)

##Inset
inset_ax = axs[ii,jj].inset_axes([0.18, 0.14, 0.4, 0.35]) 
strain_names = ['DGY1976']

for j in range(len(strain_names)):
    df2 = df_03.loc[(df_03['Strain'] == strain_names[j])]
    data = df2[['Generation', 'Multicopy_mCitrine']].to_numpy()
    inset_ax.plot(data[:, 0], data[:, 1], '.-', lw=2, markersize=0.5)
    inset_ax.set_facecolor('khaki')
    
#Mark the generation at which the strain was isolated
    gen_iso = df2.loc[df2['Generation'] == 62]
    inset_ax.scatter(gen_iso['Generation'], gen_iso['Multicopy_mCitrine'], marker='o', facecolor = 'k', edgecolor='k', s=10, zorder=10)
    
inset_ax.set_ylim(-5, 110)
inset_ax.set_yticks([0,100], [0,100], fontsize=6.5)
inset_ax.set_xticks([0,250], [0,250], fontsize=6.5)

#-----------

ii=1
jj=0

strain_names = ['DGY1976_1','DGY1976_2','DGY1976_3','DGY1976_4']

shift_values = np.linspace(-3.5, 3.5, len(strain_names))  

for j in range(len(strain_names)):
    df2 = df_2.loc[df_2['Strain'] == strain_names[j]]
    data = df2[['Generation','Multicopy_mCitrine']].to_numpy()
    axs[ii, jj].plot(data[:, 0], data[:, 1] + shift_values[j], '.-', color=gb_cols[j], lw=2)
    
axs[ii,jj].text(80,115,'G_8')
axs[ii,jj].set_ylim(-5,110)

##Inset
inset_ax = axs[ii,jj].inset_axes([0.18, 0.14, 0.4, 0.35]) 
strain_names = ['DGY1976']

for j in range(len(strain_names)):
    df2 = df_03.loc[(df_03['Strain'] == strain_names[j])]
    data = df2[['Generation', 'Multicopy_mCitrine']].to_numpy()
    inset_ax.plot(data[:, 0], data[:, 1], '.-', lw=2, markersize=0.5)
    inset_ax.set_facecolor('khaki')
    
#Mark the generation at which the strain was isolated
    gen_iso = df2.loc[df2['Generation'] == 62]
    inset_ax.scatter(gen_iso['Generation'], gen_iso['Multicopy_mCitrine'], marker='o', facecolor = 'k', edgecolor='k', s=10, zorder=10)
    
inset_ax.set_ylim(-5, 110)
inset_ax.set_yticks([0,100], [0,100], fontsize=6.5)
inset_ax.set_xticks([0,250], [0,250], fontsize=6.5)

#-----------

ii=1
jj=1

strain_names = ['DGY1920_1','DGY1920_2','DGY1920_3','DGY1920_4']

for j in range(len(strain_names)):
    df2 = df_2.loc[df_2['Strain'] == strain_names[j]]
    data = df2[['Generation','Multicopy_mCitrine']].to_numpy()
    axs[ii, jj].plot(data[:, 0], data[:, 1], '.-', color=gb_cols[j], lw=2)
    
axs[ii,jj].text(80,115,'G_9')
axs[ii,jj].set_ylim(-5,110)

##Inset
inset_ax = axs[ii,jj].inset_axes([0.18, 0.14, 0.4, 0.35]) 
strain_names = ['DGY1883']

for j in range(len(strain_names)):
    df2 = df_03.loc[(df_03['Strain'] == strain_names[j]) & (df_03['Generation'] <= 150)]
    data = df2[['Generation', 'Multicopy_mCitrine']].to_numpy()
    inset_ax.plot(data[:, 0], data[:, 1], '.-', lw=2, markersize=0.5)
    inset_ax.set_facecolor('khaki')
    
#Mark the generation at which the strain was isolated
    gen_iso = df2.loc[df2['Generation'] == 70]
    inset_ax.scatter(gen_iso['Generation'], gen_iso['Multicopy_mCitrine'], marker='o', facecolor = 'k', edgecolor='k', s=10, zorder=10)
    
inset_ax.set_ylim(-5, 110)
inset_ax.set_xlim(-5, 160)
inset_ax.set_yticks([0,100], [0,100], fontsize=6.5)
inset_ax.set_xticks([0,150], [0,150], fontsize=6.5)

#-----------

ii=1
jj=2

strain_names = ['DGY2702_1','DGY2702_2','DGY2702_3','DGY2702_4']

for j in range(len(strain_names)):
    df2 = df_2.loc[df_2['Strain'] == strain_names[j]]
    data = df2[['Generation','Multicopy_mCherry','1copy_mCherry']].to_numpy()
    Multicopy_mCherry_adj = (data[:,1] / (data[:,1] + data[:,2])) * 100
    axs[ii,jj].plot( data[:,0], Multicopy_mCherry_adj,'.-',color=r_cols[j],lw=2)
    
axs[ii,jj].text(80,115,'M_10')
axs[ii,jj].set_ylim(-5,110)

##Inset
inset_ax = axs[ii,jj].inset_axes([0.18, 0.14, 0.4, 0.35])
strain_names = ['DGY2702']

for j in range(len(strain_names)):
    df2 = df_02.loc[df_02['Strain'] == strain_names[j]]
    data = df2[['Generation', 'Multicopy_mCherry']].to_numpy()
    inset_ax.plot(data[:, 0], data[:, 1], '.-', lw=2, markersize=0.5)
    inset_ax.set_facecolor('khaki')

#Mark the generation at which the strain was isolated
    gen_iso = df2.loc[df2['Generation'] == 265.3]
    inset_ax.scatter(gen_iso['Generation'], gen_iso['Multicopy_mCherry'], marker='o', facecolor = 'k', edgecolor='k', s=10, zorder=10)

inset_ax.set_ylim(-5, 110)
inset_ax.set_yticks([0,100], [0,100], fontsize=6.5)
inset_ax.set_xticks([0,250], [0,250], fontsize=6.5)

#-----------

ii=1
jj=3

strain_names = ['DGY2769_1','DGY2769_2','DGY2769_3','DGY2769_4']

for j in range(len(strain_names)):
    df2 = df_2.loc[df_2['Strain'] == strain_names[j]]
    data = df2[['Generation','Multicopy_mCherry']].to_numpy()
    axs[ii, jj].plot(data[:, 0], data[:, 1], '.-', color=p_cols[j], lw=2)

    
axs[ii,jj].text(80,115,'P_11')
axs[ii,jj].set_ylim(-5,110)

##Inset
inset_ax = axs[ii,jj].inset_axes([0.18, 0.14, 0.4, 0.35])
strain_names = ['DGY2769']

for j in range(len(strain_names)):
    df2 = df_02.loc[df_02['Strain'] == strain_names[j]]
    data = df2[['Generation', 'Multicopy_mCherry']].to_numpy()
    inset_ax.plot(data[:, 0], data[:, 1], '.-', lw=2, markersize=0.5)
    inset_ax.set_facecolor('khaki')

#Mark the generation at which the strain was isolated
    gen_iso = df2.loc[df2['Generation'] == 88.4]
    inset_ax.scatter(gen_iso['Generation'], gen_iso['Multicopy_mCherry'], marker='o', facecolor = 'k', edgecolor='k', s=10, zorder=10)

inset_ax.set_ylim(-5, 110)
inset_ax.set_yticks([0,100], [0,100], fontsize=6.5)
inset_ax.set_xticks([0,250], [0,250], fontsize=6.5)

#-----------

axs[ii,jj].set_xlim(-5,225)   
axs[1,0].set_ylabel('                                       Percentage of cells with CNV')
axs[1,1].set_xlabel('                                    Generations')

plt.tight_layout()

#plt.savefig('Fig_3B.png', pad_inches= 0.05, dpi=300)
#plt.savefig('Fig_3B.pdf', pad_inches= 0.05, dpi=300)

plt.show()

#-----------------------------------

###Third panel

fig, axs = plt.subplots(2,4, sharex=True, sharey=True, figsize=(8,4))
fig.subplots_adjust(wspace=0.1,hspace=0.)

gb_cols = ['limegreen', 'green', 'darkturquoise', 'royalblue']
r_cols = ['peru', 'brown', 'orange', 'orangered']

#-----------

ii=0
jj=0

strain_names = ['DGY2749_1','DGY2749_2','DGY2749_3','DGY2749_4']

for j in range(len(strain_names)):
    df2 = df_2.loc[df_2['Strain'] == strain_names[j]]
    data = df2[['Generation','Multicopy_mCherry','1copy_mCherry']].to_numpy()
    Multicopy_mCherry_adj = (data[:,1] / (data[:,1] + data[:,2])) * 100
    axs[ii,jj].plot( data[:,0], Multicopy_mCherry_adj,'.-',color=r_cols[j],lw=2)
    
axs[ii,jj].text(80,115,'GM_12')
axs[ii,jj].set_ylim(-5,110)

##Inset
inset_ax = axs[ii,jj].inset_axes([0.55, 0.45, 0.4, 0.35]) 

strain_names = ['DGY2748']

for j in range(len(strain_names)):
    df2 = df_02.loc[df_02['Strain'] == strain_names[j]]
    data = df2[['Generation', 'Multicopy_mCherry']].to_numpy()
    inset_ax.plot(data[:, 0], data[:, 1], '.-', lw=2, markersize=0.5, zorder=0)
    inset_ax.set_facecolor('khaki')
    
#Mark the generation at which the strain was isolated
    gen_iso = df2.loc[df2['Generation'] == 176.8]
    inset_ax.scatter(gen_iso['Generation'], gen_iso['Multicopy_mCherry'], marker='o', facecolor = 'k', edgecolor='k', s=10, zorder=10)
    
inset_ax.set_ylim(-5, 110)
inset_ax.set_xticks([0,250], [0,250], fontsize=6.5)
inset_ax.set_yticks([0,100], [0,100], fontsize=6.5)
        
#-----------

ii=1
jj=0

strain_names = ['DGY2749_1','DGY2749_2','DGY2749_3','DGY2749_4']

for j in range(len(strain_names)):
    df2 = df_2.loc[df_2['Strain'] == strain_names[j]]
    data = df2[['Generation','Multicopy_mCitrine','1copy_mCitrine']].to_numpy()
    Multicopy_mCitrine_adj = (data[:,1] / (data[:,1] + data[:,2])) * 100
    axs[ii,jj].plot( data[:,0], Multicopy_mCitrine_adj,'.-',color=gb_cols[j],lw=2)
    
axs[ii,jj].set_ylim(-5,110)

##Inset
inset_ax = axs[ii,jj].inset_axes([0.18, 0.12, 0.4, 0.35]) 

strain_names = ['DGY2748']

for j in range(len(strain_names)):
    df2 = df_02.loc[df_02['Strain'] == strain_names[j]]
    data = df2[['Generation', 'Multicopy_mCitrine']].to_numpy()
    inset_ax.plot(data[:, 0], data[:, 1], '.-', lw=2, markersize=0.5)
    inset_ax.set_facecolor('khaki')
    
#Mark the generation at which the strain was isolated
    gen_iso = df2.loc[df2['Generation'] == 176.8]
    inset_ax.scatter(gen_iso['Generation'], gen_iso['Multicopy_mCitrine'], marker='o', facecolor = 'k', edgecolor='k', s=10, zorder=10)
    
inset_ax.set_ylim(-5, 110)
inset_ax.set_xticks([0,250], [0,250], fontsize=6.5)
inset_ax.set_yticks([0,100], [0,100], fontsize=6.5)

#-----------

ii=0
jj=1

strain_names = ['DGY2750_1','DGY2750_2','DGY2750_3','DGY2750_4']

for j in range(len(strain_names)):
    df2 = df_2.loc[df_2['Strain'] == strain_names[j]]
    data = df2[['Generation','Multicopy_mCherry','1copy_mCherry']].to_numpy()
    Multicopy_mCherry_adj = (data[:,1] / (data[:,1] + data[:,2])) * 100
    axs[ii,jj].plot( data[:,0], Multicopy_mCherry_adj,'.-',color=r_cols[j],lw=2)
    
axs[ii,jj].text(80,115,'GM_13')
axs[ii,jj].set_ylim(-5,110)

##Inset
inset_ax = axs[ii,jj].inset_axes([0.55, 0.45, 0.4, 0.35]) 

strain_names = ['DGY2748']

for j in range(len(strain_names)):
    df2 = df_02.loc[df_02['Strain'] == strain_names[j]]
    data = df2[['Generation', 'Multicopy_mCherry']].to_numpy()
    inset_ax.plot(data[:, 0], data[:, 1], '.-', lw=2, markersize=0.5)
    inset_ax.set_facecolor('khaki')
    
#Mark the generation at which the strain was isolated
    gen_iso = df2.loc[df2['Generation'] == 176.8]
    inset_ax.scatter(gen_iso['Generation'], gen_iso['Multicopy_mCherry'], marker='o', facecolor = 'k', edgecolor='k', s=10, zorder=10)
    
inset_ax.set_ylim(-5, 110)
inset_ax.set_xticks([0,250], [0,250], fontsize=6.5)
inset_ax.set_yticks([0,100], [0,100], fontsize=6.5)

#-----------

ii=1
jj=1

strain_names = ['DGY2750_1','DGY2750_2','DGY2750_3','DGY2750_4']

for j in range(len(strain_names)):
    df2 = df_2.loc[df_2['Strain'] == strain_names[j]]
    data = df2[['Generation','Multicopy_mCitrine','1copy_mCitrine']].to_numpy()
    Multicopy_mCitrine_adj = (data[:,1] / (data[:,1] + data[:,2])) * 100
    axs[ii,jj].plot( data[:,0], Multicopy_mCitrine_adj,'.-',color=gb_cols[j],lw=2)
    
axs[ii,jj].set_ylim(-5,110)

##Inset
inset_ax = axs[ii,jj].inset_axes([0.18, 0.12, 0.4, 0.35], zorder=-1) 

strain_names = ['DGY2748']

for j in range(len(strain_names)):
    df2 = df_02.loc[df_02['Strain'] == strain_names[j]]
    data = df2[['Generation', 'Multicopy_mCitrine']].to_numpy()
    inset_ax.plot(data[:, 0], data[:, 1], '.-', lw=2, markersize=0.5, zorder=0)
    inset_ax.set_facecolor('khaki')
    
#Mark the generation at which the strain was isolated
    gen_iso = df2.loc[df2['Generation'] == 176.8]
    inset_ax.scatter(gen_iso['Generation'], gen_iso['Multicopy_mCitrine'], marker='o', facecolor = 'k', edgecolor='k', s=10, zorder=10)
    
inset_ax.set_ylim(-5, 110)
inset_ax.set_xticks([0,250], [0,250], fontsize=6.5)
inset_ax.set_yticks([0,100], [0,100], fontsize=6.5)

#-----------

ii=0
jj=2

strain_names = ['DGY2751_2','DGY2751_3','DGY2751_4']

for j in range(len(strain_names)):
    df2 = df_2.loc[df_2['Strain'] == strain_names[j]]
    data = df2[['Generation','Multicopy_mCherry','1copy_mCherry']].to_numpy()
    Multicopy_mCherry_adj = (data[:,1] / (data[:,1] + data[:,2])) * 100
    axs[ii,jj].plot( data[:,0], Multicopy_mCherry_adj,'.-',color=r_cols[j],lw=2)
        
axs[ii,jj].text(80,115,'GM_14')
axs[ii,jj].set_ylim(-5,110)

##Inset
inset_ax = axs[ii,jj].inset_axes([0.55, 0.45, 0.4, 0.35]) 

strain_names = ['DGY2748']

for j in range(len(strain_names)):
    df2 = df_02.loc[df_02['Strain'] == strain_names[j]]
    data = df2[['Generation', 'Multicopy_mCherry']].to_numpy()
    inset_ax.plot(data[:, 0], data[:, 1], '.-', lw=2, markersize=0.5)
    inset_ax.set_facecolor('khaki')
    
#Mark the generation at which the strain was isolated
    gen_iso = df2.loc[df2['Generation'] == 176.8]
    inset_ax.scatter(gen_iso['Generation'], gen_iso['Multicopy_mCherry'], marker='o', facecolor = 'k', edgecolor='k', s=10, zorder=10)
    
inset_ax.set_ylim(-5, 110)
inset_ax.set_xticks([0,250], [0,250], fontsize=6.5)
inset_ax.set_yticks([0,100], [0,100], fontsize=6.5)

#-----------

ii=1
jj=2

strain_names = ['DGY2751_2','DGY2751_3','DGY2751_4']

for j in range(len(strain_names)):
    df2 = df_2.loc[df_2['Strain'] == strain_names[j]]
    data = df2[['Generation','Multicopy_mCitrine','1copy_mCitrine']].to_numpy()
    Multicopy_mCitrine_adj = (data[:,1] / (data[:,1] + data[:,2])) * 100
    axs[ii,jj].plot( data[:,0], Multicopy_mCitrine_adj,'.-',color=gb_cols[j],lw=2)
    
axs[ii,jj].set_ylim(-5,110)

##Inset
inset_ax = axs[ii,jj].inset_axes([0.18, 0.12, 0.4, 0.35], zorder=-1) 

strain_names = ['DGY2748']

for j in range(len(strain_names)):
    df2 = df_02.loc[df_02['Strain'] == strain_names[j]]
    data = df2[['Generation', 'Multicopy_mCitrine']].to_numpy()
    inset_ax.plot(data[:, 0], data[:, 1], '.-', lw=2, markersize=0.5)
    inset_ax.set_facecolor('khaki')
    
#Mark the generation at which the strain was isolated
    gen_iso = df2.loc[df2['Generation'] == 176.8]
    inset_ax.scatter(gen_iso['Generation'], gen_iso['Multicopy_mCitrine'], marker='o', facecolor = 'k', edgecolor='k', s=10, zorder=10)
    
inset_ax.set_ylim(-5, 110)
inset_ax.set_xticks([0,250], [0,250], fontsize=6.5)
inset_ax.set_yticks([0,100], [0,100], fontsize=6.5)

#-----------

ii=0
jj=3

strain_names = ['DGY2755_1','DGY2755_2','DGY2755_3','DGY2755_4']

for j in range(len(strain_names)):
    df2 = df_2.loc[df_2['Strain'] == strain_names[j]]
    data = df2[['Generation','Multicopy_mCherry','1copy_mCherry']].to_numpy()
    Multicopy_mCherry_adj = (data[:,1] / (data[:,1] + data[:,2])) * 100
    axs[ii,jj].plot( data[:,0], Multicopy_mCherry_adj,'.-',color=r_cols[j],lw=2)
    
axs[ii,jj].text(80,115,'GM_15')
axs[ii,jj].set_ylim(-5,110)

##Inset
inset_ax = axs[ii,jj].inset_axes([0.55, 0.45, 0.4, 0.35]) 

strain_names = ['DGY2748']

for j in range(len(strain_names)):
    df2 = df_02.loc[df_02['Strain'] == strain_names[j]]
    data = df2[['Generation', 'Multicopy_mCherry']].to_numpy()
    inset_ax.plot(data[:, 0], data[:, 1], '.-', lw=2, markersize=0.5)
    inset_ax.set_facecolor('khaki')
    
#Mark the generation at which the strain was isolated
    gen_iso_2 = df2.loc[df2['Generation'] == 147.4]
    inset_ax.scatter(gen_iso_2['Generation'], gen_iso_2['Multicopy_mCherry'], marker='o', facecolor = 'k', edgecolor='k', s=10, zorder=10)

inset_ax.set_ylim(-5, 110)
inset_ax.set_xticks([0,250], [0,250], fontsize=6.5)
inset_ax.set_yticks([0,100], [0,100], fontsize=6.5)

#-----------

ii=1
jj=3

strain_names = ['DGY2755_1','DGY2755_2','DGY2755_3','DGY2755_4']

for j in range(len(strain_names)):
    df2 = df_2.loc[df_2['Strain'] == strain_names[j]]
    data = df2[['Generation','Multicopy_mCitrine','1copy_mCitrine']].to_numpy()
    Multicopy_mCitrine_adj = (data[:,1] / (data[:,1] + data[:,2])) * 100
    axs[ii,jj].plot( data[:,0], Multicopy_mCitrine_adj,'.-',color=gb_cols[j],lw=2)
    
axs[ii,jj].set_ylim(-5,110)

##Inset
inset_ax = axs[ii,jj].inset_axes([0.18, 0.12, 0.4, 0.35]) 

strain_names = ['DGY2748']

for j in range(len(strain_names)):
    df2 = df_02.loc[df_02['Strain'] == strain_names[j]]
    data = df2[['Generation', 'Multicopy_mCitrine']].to_numpy()
    inset_ax.plot(data[:, 0], data[:, 1], '.-', lw=2, markersize=0.5)
    inset_ax.set_facecolor('khaki')
    
#Mark the generation at which the strain was isolated 
    gen_iso_2 = df2.loc[df2['Generation'] == 147.4]
    inset_ax.scatter(gen_iso_2['Generation'], gen_iso_2['Multicopy_mCitrine'], marker='o', facecolor = 'k', edgecolor='k', s=10, zorder=10)

inset_ax.set_ylim(-5, 110)
inset_ax.set_xticks([0,250], [0,250], fontsize=6.5)
inset_ax.set_yticks([0,100], [0,100], fontsize=6.5)

#-----------

axs[ii,jj].set_xlim(-5,225)
axs[1,0].set_ylabel('                                   Percentage of cells with CNV')
axs[1,1].set_xlabel('                                    Generations')

plt.tight_layout(h_pad=-0.4) #reduces spacing between top and bottom rows

#plt.savefig('Fig_3C.png', pad_inches= 0.05, dpi=300)
#plt.savefig('Fig_3C.pdf', pad_inches= 0.05, dpi=300)

plt.show()

#-----------------------------------

