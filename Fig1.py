#Load libraries

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.special import erf

#-----------------------------------
#-----------------------------------
#-----------------------------------

###Fig 1A: Schematic of the main question 

fig = plt.figure(figsize=(5,4))

gs = fig.add_gridspec(4,2, wspace=0.25)
ax1 = fig.add_subplot(gs[1:3, 0])
ax2 = fig.add_subplot(gs[0:2, 1])
ax3 = fig.add_subplot(gs[2:4, 1])


x = np.linspace(0,1,60)
y = 0.5 * (1. + erf((x-0.4)/0.1) )            + np.random.normal(loc=0.0, scale=0.02, size=len(x))
ax1.plot(x,y,lw=2)

y2 = x*0. + 1.                                + np.random.normal(loc=0.0, scale=0.02, size=len(x))
ax2.plot(x, y2,lw=2)

y3 = 1. -  0.5 * (1. + erf((x-0.55)/0.15) )   + np.random.normal(loc=0.0, scale=0.02, size=len(x))
ax3.plot(x,y3,lw=2)

for ax in [ax1,ax2,ax3]:
    ax.set_ylim(-0.1,1.1)
    ax.set_xlim(0,1)
    ax.set_xticks([])
    ax.set_yticks([0,0.5,1],['0','50','100'])

ax1.set_ylabel('Percentage of cells with CNV')

for ax in [ax1,ax3]:
    ax.set_xlabel('Generations')
    
    
ax3.annotate('', xy=(0.48, 0.53), xytext=(0.56, 0.6), xycoords='figure fraction', textcoords='figure fraction',
                  arrowprops=dict(arrowstyle='<-,  head_width=0.3', color='k', linewidth=1.2, linestyle='-'))

ax3.annotate('', xy=(0.48, 0.38), xytext=(0.56, 0.3), xycoords='figure fraction', textcoords='figure fraction',
                  arrowprops=dict(arrowstyle='<-,  head_width=0.3', color='k', linewidth=1.2, linestyle='-'))

ax1.text(0.45,0.2,'Selection\nfor CNVs')
ax2.text(0.25,0.2,'Stable CNVs')
ax3.text(0.25,0.2,'Loss\nof CNVs')

#ax1.set_title('Selection pressure')
#ax2.set_title('Selection relaxed')

ax1.fill_between([0.15,1.1],[1.1,1.1],[-0.1,-0.1],color='khaki',zorder=0)

for ax in [ax2,ax3]:
    ax.fill_between([0.,0.2],[1.1,1.1],[-0.1,-0.1],color='khaki',zorder=0)
    
plt.tight_layout()

#plt.savefig('Fig_1A.png',dpi=300,pad_inches=0.05)
#plt.savefig('Fig_1A.pdf',dpi=300,pad_inches=0.05)

plt.show()

#-----------------------------------
#-----------------------------------
#-----------------------------------

###Fig 1B: Schematic of all CNV strains in the study 

df = pd.read_csv('All_strain_schematic_GAP1_20250122.CSV')

gridspec = dict(hspace=0.3, height_ratios=[12,0.3,3,0.3,3,1,4])
fig, axs = plt.subplots(nrows=7, ncols=1, figsize=(8,10), gridspec_kw=gridspec)
axs[1].set_visible(False)
axs[3].set_visible(False)
axs[5].set_visible(False)

ax1 = axs[0]
ax2 = axs[2]
ax3 = axs[4]
ax4 = axs[6]

#print(np.max(df['CNV_Endpoint']/1000))


color_array = ['darkgreen','springgreen','springgreen','springgreen','springgreen','springgreen','springgreen','greenyellow','greenyellow']
pattern = ['xxx','///','///','///','///','///','///','','']
for i in range(len(df['Strain'])):
    ax1.fill_between([df['CNV_Startpoint'][i]/1000,df['CNV_Endpoint'][i]/1000],[-i-0.3,-i-0.3],[-i+0.3,-i+0.3],  color=color_array[i], hatch=pattern[i], edgecolor="k", linewidth=1.0)
    #ax.plot([df['CNV_Startpoint'][i]/1000,df['CNV_Endpoint'][i]/1000],[i,i], ms=10, color=color_array[i])
#plt.yscale('log')
ax1.set_xticks(ticks=[400,500,600,666],labels=['0','500','600','666'], color='green');
ax1.set_yticks(ticks=-1. * np.arange(0,9,1),labels=df['Names']);
ax1.set_xlabel('Chromosome XI coordinate (kb)',fontsize=12, color='green')
#ax1.set_ylabel('Chr 11',fontsize=14)
ax1.set_ylim(-len(df['Strain'])-0.3,1.3)
#ax1.set_xlim(400, np.max(df['CNV_Endpoint']/1000))
ax1.set_xlim(400, 666)

ax1.plot([516,516],[-100,100],':k',lw=2)
ax1.tick_params(axis='x', which='major', labelsize=14, color='green')


#----- legend box, and text
color_array = ['darkgreen','springgreen','greenyellow','indianred','pink']
pattern = ['xxx','///','','///','']

label_name = ['4-copy GAP1','3-copy GAP1','2-copy GAP1']

for i in range(3):
    ax1.fill_between([585,595],[-i/1.2-6-0.2,-i/1.2-6-0.2],[-i/1.2+0.2-6,-i/1.2+0.2-6],  color=color_array[i], hatch=pattern[i], edgecolor="k", linewidth=1.0)
    #ax.plot([df['CNV_Startpoint'][i]/1000,df['CNV_Endpoint'][i]/1000],[i,i], ms=10, color=color_array[i])
    ax1.text(600, -i/1.2-0.2 -6, label_name[i],fontsize=14)

ax1.plot([580,660],[-5.25,-5.25],'-k')
ax1.plot([580,660],[-8.25,-8.25],'-k')
ax1.plot([580,580],[-5.25,-8.25],'-k')
ax1.plot([660,660],[-5.25,-8.25],'-k')
#--------------------------------------



df = pd.read_csv('All_strain_schematic_MEP2.CSV')
color_array = ['orangered']
pattern = ['///']
for i in range(len(df['Strain'])):
    ax2.fill_between([df['CNV_Startpoint'][i]/1000,df['CNV_Endpoint'][i]/1000],[-i-0.25,-i-0.25],[-i+0.25,-i+0.25],  color=color_array[i], hatch=pattern[i], edgecolor="k", linewidth=1.0)
    #ax.plot([df['CNV_Startpoint'][i]/1000,df['CNV_Endpoint'][i]/1000],[i,i], ms=10, color=color_array[i])
#plt.yscale('log')
ax2.set_xticks(ticks=[230, 300, 400, 496],labels=['0','300','400','784'], color='red');
ax2.set_yticks(ticks=-1. * np.arange(0,len(df['Strain']),1),labels=df['Names']);
ax2.set_xlabel('Chromosome XIV coordinate (kb)',fontsize=12, color='red')
ax2.set_ylabel('      Strain name',fontsize=17)
ax2.set_ylim(-len(df['Strain'])-0.3,1.3)
ax2.set_xlim(230, 496) #np.max(df['CNV_Endpoint']/1000))

ax2.plot([359,359],[-100,100],':k',lw=2)
ax2.tick_params(axis='x', which='major', labelsize=14, color='red')

#----- legend box, and text
color_array = ['orangered','orange']
pattern = ['///','']

label_name = ['3-copy MEP2','2-copy MEP2']

for i in range(2):
    ax2.fill_between([415,425],[-i/1.2-0.2 + 0.5/1.2, -i/1.2-0.2 + 0.5/1.2],[-i/1.2+0.2 + 0.5/1.2, -i/1.2+0.2 + 0.5/1.2],  color=color_array[i], hatch=pattern[i], edgecolor="k", linewidth=1.0)
    #ax.plot([df['CNV_Startpoint'][i]/1000,df['CNV_Endpoint'][i]/1000],[i,i], ms=10, color=color_array[i])
    ax2.text(430, -i/1.2-0.2 + 0.5/1.2, label_name[i],fontsize=14)

ax2.plot([410,490],[1,1],'-k')
ax2.plot([410,490],[-1,-1],'-k')
ax2.plot([410,410],[1,-1],'-k')
ax2.plot([490,490],[1,-1],'-k')
#--------------------------------------



df = pd.read_csv('All_strain_schematic_PUT4.CSV')
color_array = ['violet']
pattern = ['///']

label_name = ['3-copy PUT4']

for i in range(len(df['Strain'])):
    ax3.fill_between([df['CNV_Startpoint'][i]/1000,df['CNV_Endpoint'][i]/1000],[-i-0.25,-i-0.25],[-i+0.25,-i+0.25],  color=color_array[i], hatch=pattern[i], edgecolor="k", linewidth=1.0)
    #ax.plot([df['CNV_Startpoint'][i]/1000,df['CNV_Endpoint'][i]/1000],[i,i], ms=10, color=color_array[i])
#plt.yscale('log')
ax3.set_xticks(ticks=[850,925,1025,1116],labels=['0','925','1025','1091']);
ax3.set_yticks(ticks=-1. * np.arange(0,len(df['Strain']),1),labels=df['Names']);
ax3.set_xlabel('Chromosome XV coordinate (kb)',fontsize=12)
#ax3.set_ylabel('Chr 15',fontsize=14)
ax3.set_ylim(-len(df['Strain'])-0.3,1.3)
ax3.set_xlim(850, 1116) #np.max(df['CNV_Endpoint']/1000))

ax3.plot([987,987],[-100,100],':k',lw=2)

#----- legend box, and text
color_array = ['violet']
pattern = ['///']

label_name = ['3-copy PUT4']

for i in range(1):
    ax3.fill_between([1035,1045],[-i/1.2-0.2, -i/1.2-0.2],[-i/1.2+0.2, -i/1.2+0.2],  color=color_array[i], hatch=pattern[i], edgecolor="k", linewidth=1.0)
    #ax.plot([df['CNV_Startpoint'][i]/1000,df['CNV_Endpoint'][i]/1000],[i,i], ms=10, color=color_array[i])
    ax3.text(1050, -i/1.2-0.2, label_name[i],fontsize=14)

ax3.plot([1030,1110],[0.5,0.5],'-k')
ax3.plot([1030,1110],[-0.5,-0.5],'-k')
ax3.plot([1030,1030],[0.5,-0.5],'-k')
ax3.plot([1110,1110],[0.5,-0.5],'-k')
#--------------------------------------

for i in range(1):
    color ='orange'
    ax4.fill_between([0,1],[-i+0.3,-i+0.3],[-i+0.05,-i+0.05],  color=color, edgecolor="k", linewidth=1.0)
    color ='greenyellow'
    ax4.fill_between([0,1],[-i-0.3,-i-0.3],[-i-0.05,-i-0.05],  color=color, edgecolor="k", linewidth=1.0)

ax4.set_xticks(ticks=[0,0.3333,0.6666,1],labels=['0','222','444','666'], color='green');
ax4.set_yticks(ticks=-1. * np.arange(0,1,1),labels=['GM_{12 to 15}'],rotation=90,verticalalignment='center');
ax4.set_xlabel('Chromosome      or        coordinate (kb)',fontsize=12)
#ax4.set_ylabel('Chr 11, 14\n',fontsize=14)
ax4.set_ylim(-1,1)
ax4.set_xlim(0, 1) #np.max(df['CNV_Endpoint']/1000))

ax4.plot([359/784,359/784],[0,100],':k',lw=2)
ax4.plot([516/666,516/666],[-100,0],':k',lw=2)

for ax in axs[:6]:
    ax.tick_params(axis='both', which='major', labelsize=14,
                    direction='in',length=8, width=1.3, bottom=True, left=True, top=False, right=False, zorder=15.0)
    #ax.tick_params(axis='both', which='minor', labelsize=16,
                    #direction='in',length=4, width=1.3, bottom=True, left=True, top=True, right=True, zorder=15.0)
    ax.spines['left'].set_linewidth(1.3)
    ax.spines['right'].set_linewidth(1.3)
    ax.spines['top'].set_linewidth(1.3)
    ax.spines['bottom'].set_linewidth(1.3)
    
ax4.tick_params(axis='x', which='major', labelsize=14, color='green',
                direction='in',length=8, width=1.3, bottom=True, left=True, top=False, right=False, zorder=15.0)
ax4.tick_params(axis='y', which='major', labelsize=14,
                direction='in',length=0, width=0, bottom=True, left=True, top=False, right=False, zorder=15.0)
#ax.tick_params(axis='both', which='minor', labelsize=16,
                #direction='in',length=4, width=1.3, bottom=True, left=True, top=True, right=True, zorder=15.0)

ax4.spines['left'].set_linewidth(1.3)
ax4.spines['right'].set_linewidth(1.3)
#ax4.spines['top'].set_linewidth(1.3)
#ax4.spines['bottom'].set_linewidth(1.3)


ax4_twin = ax4.twiny()
ax4_twin.set_xticks(ticks=[0,0.3333,0.6666,1],labels=['0','262','523','784'],color='red');
ax4_twin.tick_params(axis='x', which='major', labelsize=14, color='red',
                direction='in',length=8, width=1.3, bottom=False, left=False, top=True, right=False, zorder=15.0)


ax4.spines['bottom'].set_color('green')
ax4_twin.spines['top'].set_color('red')
ax4.spines['top'].set_lw(1.3)
ax4.spines['bottom'].set_lw(1.3)

ax4.text(0.43,-1.8,'XI',color='green',fontsize=12)
ax4.text(0.5,-1.8,'XIV',color='red',fontsize=12)

#plt.savefig('Fig_1B.png',bbox_inches='tight',pad_inches = 0.05,dpi=300)
#plt.savefig('Fig_1B.pdf',bbox_inches='tight',pad_inches = 0.05,dpi=300)

plt.show()

#-----------------------------------
#-----------------------------------
#-----------------------------------
