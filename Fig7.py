###Load libraries

import pandas as pd
import pysam
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import median_filter, gaussian_filter1d, label
import FlowCal
from matplotlib.gridspec import GridSpec

#------------------------------------------------------

###Read in bedgraph files

bedgraph_1 = "PacBio_1886.bedgraph"
bedgraph_1886 = pd.read_csv(bedgraph_1, sep="\t", header=None, names=["chrom", "start", "end", "value"])

bedgraph_2 = "PacBio_3417.bedgraph"
bedgraph_3417 = pd.read_csv(bedgraph_2, sep="\t", header=None, names=["chrom", "start", "end", "value"])

bedgraph_3 = "PacBio_3410.bedgraph"
bedgraph_3410 = pd.read_csv(bedgraph_3, sep="\t", header=None, names=["chrom", "start", "end", "value"])

bedgraph_4 = "PacBio_3414.bedgraph"
bedgraph_3414 = pd.read_csv(bedgraph_4, sep="\t", header=None, names=["chrom", "start", "end", "value"])

bedgraph_5 = "PacBio_3418.bedgraph"
bedgraph_3418 = pd.read_csv(bedgraph_5, sep="\t", header=None, names=["chrom", "start", "end", "value"])

bedgraph_6 = "PacBio_3420.bedgraph"
bedgraph_3420 = pd.read_csv(bedgraph_6, sep="\t", header=None, names=["chrom", "start", "end", "value"])

#------------------------------------------------------

###Figure 7: Plot read depth and fluorescence of all strains

fig = plt.figure(figsize=(7,4),layout="constrained")

gs = GridSpec(3, 2, figure=fig)
ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[1,0])
ax3 = fig.add_subplot(gs[2,0])

ax4 = fig.add_subplot(gs[:,1])

chrom = "NC_001143.9"
start = 400000
end = 650000

strains = [
    ("CNV strain", bedgraph_1886),
    ("Partial\nrevertant", bedgraph_3417), 
    ("Comp", bedgraph_3410),
    ("         lete", bedgraph_3414),
    ("        \nrever", bedgraph_3418),
    ("         \n        tants", bedgraph_3420)
]

colors = ["darkgreen", "lime", "tab:cyan", "tab:orange", "tab:brown", "tab:olive"]
linestyles = ["-", "-", "-", "--", "-.", ":"]

for i, (title, bedgraph) in enumerate(strains):
    chrom_subset = bedgraph[bedgraph["chrom"] == chrom].copy()  # avoid modifying original
    
    # Compute median of the 5,000 bp region before 510,000
    preceding_median = chrom_subset.loc[
        (chrom_subset["start"] >= 500050) & (chrom_subset["start"] < 510000), "value"
    ].median()

    # Replace values between 510,000 and 525,000 with that median
    chrom_subset.loc[
        (chrom_subset["start"] >= 515000) & (chrom_subset["start"] <= 518000), "value"
    ] = preceding_median


    if i==0:
        ax = ax1
    elif i==1:
        ax = ax2
    elif i==2:
        ax = ax3
    elif i==3:
        ax = ax3
    elif i==4:
        ax = ax3
    elif i==5:
        ax = ax3
    else:
        ax = ax4
   
    mask = (chrom_subset["start"] > 4.8e5) & (chrom_subset["start"] < 6.5e5)
    base_value = np.mean(chrom_subset["value"][~mask])

    sigma = 20
    lw = 3

    y = gaussian_filter1d(chrom_subset["value"] / base_value, sigma)
    ax.plot(chrom_subset["start"], y, lw=lw, color=colors[i], linestyle=linestyles[i], alpha=0.8)

    ###Split reads for CNV strain
    # Load SAM file and extract split-read positions
    split_df = pd.read_csv(
        "Split_reads_1886.sam",
        sep="\t",
        comment="@",
        header=None,
        usecols=[2,3,5],
        names=["chrom", "pos", "cigar"]
    )
    # Keep only reads on this chromosome and with split CIGAR (S or H)
    split_df = split_df[
        (split_df["chrom"] == chrom) &
        (split_df["cigar"].str.contains("S|H"))
    ]
    # Filter to the two desired regions
    split_sub = split_df[
        ((split_df["pos"] >= 497500) & (split_df["pos"] <= 500000)) |
        ((split_df["pos"] >= 634000) & (split_df["pos"] <= 635700))
    ]
    # Plot them
    ax1.vlines(split_sub["pos"], ymin=3.5, ymax=5,
              linewidth=1, color="red", alpha=0.01)
    

    ###Split reads for partial revertant
    # Load SAM file and extract split-read positions
    split_df = pd.read_csv(
        "Split_reads_3417.sam",
        sep="\t",
        comment="@",
        header=None,
        usecols=[2,3,5],
        names=["chrom", "pos", "cigar"]
    )
    # Keep only reads on this chromosome and with split CIGAR (S or H)
    split_df = split_df[
        (split_df["chrom"] == chrom) &
        (split_df["cigar"].str.contains("S|H"))
    ]
    # Filter to the two desired regions
    split_sub = split_df[
        ((split_df["pos"] >= 497500) & (split_df["pos"] <= 500000)) |
        ((split_df["pos"] >= 634000) & (split_df["pos"] <= 635700))
    ]
    # Plot them
    ax2.vlines(split_sub["pos"], ymin=3.5, ymax=5,
              linewidth=1, color="red", alpha=0.01)


   
    if i in [2, 5]:
        ax.set_xlabel("Genomic position on Chr XI (kb)")
    if i == 1:
        ax.set_ylabel("Ploidy-normalized copy number")
    ax.set_xticks([400000, 450000, 500000, 550000, 600000, 650000])
    ax.set_xticklabels([400, 450, 500, 550, 600, 650])
    if i in [0,1,6]:
        ax.set_xticklabels([])
    if i < 3:
        ax.plot([518400, 518400], [-100, 100], "--m", zorder= -2, lw=1)
        ax.plot([515500, 515500], [-100, 100], "--b", zorder= -2, lw=1)

    ax.text(4.1e5, 5.7, f"{title}", verticalalignment="top", color=colors[i])
    ax.plot([497531, 497531], [-100, 100], ":k")
    ax.plot([635729, 635729], [-100, 100], ":k")
    ax.plot([0, 800000], [1, 1], "--", lw=1,color="grey")
    ax.plot([0, 800000], [2, 2], "--", lw=1,color="grey")
    ax.plot([0, 800000], [3, 3], "--", lw=1,color="grey")
    ax.set_xlim((400000, 650000))
    ax.set_ylim(0, 6)
    ax.set_yticks([1,3],[1,3])

    if i==0:
        ax.annotate(
            text='',
            xy=(515500, 4),    # The point to annotate (data coordinates)
            xytext=(550000, 4),          # The position of the text
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2") # Arrow properties
        )
        ax.text(550000, 4, 'GAP1', verticalalignment='center', color='m')

    if i==0:
        ax.annotate(
            text='',
            xy=(513000, 5),    # The point to annotate (data coordinates)
            xytext=(545000, 5),          # The position of the text
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2") # Arrow properties
        )
        ax.text(545000, 5, 'mCitrine', verticalalignment='center', color='b')

#-----------------------------------------

ax = ax4

data = FlowCal.io.FCSData('DGY1657_1.fcs')
channels = data.channels
df = pd.DataFrame(data, columns=channels)
mask = ( df['FSC-A']<4.1e6 ) & ( df['B2-A']>1e4 )
ax.scatter(df[mask]['FSC-A'], df[mask]['B2-A'], s=1, color='black', alpha=0.5)


data = FlowCal.io.FCSData('DGY1886_2.fcs')
channels = data.channels
df = pd.DataFrame(data, columns=channels)
mask = ( df['FSC-A']<4.1e6 ) & ( df['B2-A']>1e4 )
ax.scatter(df[mask]['FSC-A'], df[mask]['B2-A'], s=2, color='darkgreen', alpha=0.5)
data = FlowCal.io.FCSData('DGY1886_3.fcs')
channels = data.channels
df = pd.DataFrame(data, columns=channels)
mask = ( df['FSC-A']<4.1e6 ) & ( df['B2-A']>1e4 )
ax.scatter(df[mask]['FSC-A'], df[mask]['B2-A'], s=2, color='darkgreen', alpha=0.5)

data = FlowCal.io.FCSData('DGY3417.fcs')
channels = data.channels
df = pd.DataFrame(data, columns=channels)
mask = ( df['FSC-A']<4.1e6 ) & ( df['B2-A']>1e4 )
ax.scatter(df[mask]['FSC-A'], df[mask]['B2-A'], s=1, color='lime', alpha=0.5)

data = FlowCal.io.FCSData('DGY3410.fcs')
channels = data.channels
df = pd.DataFrame(data, columns=channels)
mask = ( df['FSC-A']<4.1e6 ) & ( df['B2-A']>1e4 )
ax.scatter(df[mask]['FSC-A'], df[mask]['B2-A'], s=1, color='tab:cyan', alpha=0.8)

data = FlowCal.io.FCSData('DGY3414.fcs')
channels = data.channels
df = pd.DataFrame(data, columns=channels)
mask = ( df['FSC-A']<4.1e6 ) & ( df['B2-A']>1e4 )
ax.scatter(df[mask]['FSC-A'], df[mask]['B2-A'], s=1, color='tab:orange', alpha=0.6)

data = FlowCal.io.FCSData('DGY3418.fcs')
channels = data.channels
df = pd.DataFrame(data, columns=channels)
mask = ( df['FSC-A']<4.1e6 ) & ( df['B2-A']>1e4 )
ax.scatter(df[mask]['FSC-A'], df[mask]['B2-A'], s=1, color='tab:brown', alpha=0.4)

data = FlowCal.io.FCSData('DGY3420.fcs')
channels = data.channels
df = pd.DataFrame(data, columns=channels)
mask = ( df['FSC-A']<4.1e6 ) & ( df['B2-A']>1e4 )
ax.scatter(df[mask]['FSC-A'], df[mask]['B2-A'], s=1, color='tab:olive', alpha=0.2)


ax.set_xlabel('Forward scatter area')
ax.set_ylabel('log(mCitrine fluorescence)')
ax.set_yscale('log')
ax.set_ylim(2e4,2e6)
ax.set_xlim(1e4,4.2e6)

plt.gcf().align_labels()
plt.tight_layout(h_pad=-0.2)

#plt.savefig('Fig8_A.png', dpi=500)
#plt.savefig('Fig8_A_revised.pdf', dpi=900)

plt.show()

#-------------------------------------------------------
