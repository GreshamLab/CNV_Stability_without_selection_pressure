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

bedgraph_1 = "1706.bedgraph"
bedgraph_1706 = pd.read_csv(bedgraph_1, sep="\t", header=None, names=["chrom", "start", "end", "value"])

bedgraph_2 = "1886.bedgraph"
bedgraph_1886 = pd.read_csv(bedgraph_2, sep="\t", header=None, names=["chrom", "start", "end", "value"])

bedgraph_3 = "3417.bedgraph"
bedgraph_3417 = pd.read_csv(bedgraph_3, sep="\t", header=None, names=["chrom", "start", "end", "value"])

bedgraph_4 = "3286.bedgraph"
bedgraph_3286 = pd.read_csv(bedgraph_4, sep="\t", header=None, names=["chrom", "start", "end", "value"])

bedgraph_5 = "3410.bedgraph"
bedgraph_3410 = pd.read_csv(bedgraph_5, sep="\t", header=None, names=["chrom", "start", "end", "value"])

#------------------------------------------------------

###Input bam files and compute coverage

def compute_coverage(bam_file, chrom, start, end):
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        return np.array([bam.count(chrom, pos, pos+1) for pos in range(start, end)])

chrom = "NC_001143.9"
start = 400000
end = 650000

coverage_1706 = compute_coverage("1706_split.bam", chrom, start, end)
coverage_1886 = compute_coverage("1886_split.bam", chrom, start, end)
coverage_3417 = compute_coverage("3417_split.bam", chrom, start, end)
coverage_3286 = compute_coverage("3286_split.bam", chrom, start, end)
coverage_3410 = compute_coverage("3410_split.bam", chrom, start, end)

#------------------------------------------------------

###Figure 7: Plot read depth and fluorescence of all strains

# Set coverage to zero for the specified region
coverage_1706[(np.arange(start, end) >= 510000) & (np.arange(start, end) <= 525000)] = 0
coverage_1886[(np.arange(start, end) >= 510000) & (np.arange(start, end) <= 525000)] = 0
coverage_3417[(np.arange(start, end) >= 510000) & (np.arange(start, end) <= 525000)] = 0
coverage_3286[(np.arange(start, end) >= 510000) & (np.arange(start, end) <= 525000)] = 0
coverage_3410[(np.arange(start, end) >= 510000) & (np.arange(start, end) <= 525000)] = 0

fig = plt.figure(figsize=(8,6),layout="constrained")

gs = GridSpec(5, 2, figure=fig)
ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[1,0])
ax3 = fig.add_subplot(gs[2,0])
ax4 = fig.add_subplot(gs[3,0])
ax5 = fig.add_subplot(gs[4,0])

ax6 = fig.add_subplot(gs[:,1])

strains = [
    ("Ancestor", bedgraph_1706, coverage_1706),
    ("CNV strain", bedgraph_1886, coverage_1886),
    ("Partial\nrevertant", bedgraph_3417, coverage_3417),
    ("Complete\nrevertant A", bedgraph_3286, coverage_3286),
    ("Complete\nrevertant B", bedgraph_3410, coverage_3410),
]

#Bedgraph plots

window_size = 5000  # Rolling median window in bp

colors = ["tab:blue", "darkgreen", "lime", "tab:olive", "tab:brown"]

for i, (title, bedgraph, coverage) in enumerate(strains):
    chrom_subset = bedgraph[bedgraph["chrom"] == chrom]

    if i==0:
        ax = ax1
    elif i==1:
        ax = ax2
    elif i==2:
        ax = ax3
    elif i==3:
        ax = ax4
    else:
        ax = ax5
   
    mask = (chrom_subset["start"] > 4.8e5) & (chrom_subset["start"] < 6.5e5)
    base_value = np.mean(chrom_subset["value"][~mask])

    sigma = 50
    lw = 1 if i == 0 else 0.7

    # Compute rolling median
    rolling_median = median_filter(chrom_subset["value"] / base_value, size=window_size)

    ax.plot(chrom_subset["start"], gaussian_filter1d(chrom_subset["value"] / base_value, sigma), lw=lw, color=colors[i], alpha=0.8)
    ax.plot(chrom_subset["start"], rolling_median, lw=1, color="black") 


   
    if i == 4:
        ax.set_xlabel("Genomic position on Chr 11 (bp)")
        ax.set_xticks([400000, 450000, 500000, 550000, 600000, 650000]) 
    if i == 2:
        ax.set_ylabel("Normalised read depth")
    if i == 0:
        ax.set_title("Copy number from sequencing")
    if i != 4:
        ax.set_xticklabels([])

   
    ax2.vlines([497531, 635729], ymin=0, ymax=3, colors='red', linewidth=2)
    ax3.vlines([497531, 635729], ymin=0, ymax=2, colors='red', linewidth=2)

    ax.text(4.1e5, 5.5, f"{title}", verticalalignment="top", color=colors[i])
    ax.plot([497531, 497531], [-100, 100], ":k")
    ax.plot([635729, 635729], [-100, 100], ":k")
    ax.plot([0, 800000], [1, 1], "--", lw=1,color="grey")
    ax.plot([0, 800000], [2, 2], "--", lw=1,color="grey")
    ax.plot([0, 800000], [3, 3], "--", lw=1,color="grey")
    ax.set_xlim((400000, 650000))
    ax.set_ylim(0, 6)
    ax.set_yticks([1,3],[1,3])




ax = ax6

data = FlowCal.io.FCSData('DGY1657_1.fcs')
channels = data.channels
df = pd.DataFrame(data, columns=channels)
mask = ( df['FSC-A']<4.1e6 ) & ( df['B2-A']>1e4 )
ax.scatter(df[mask]['FSC-A'], df[mask]['B2-A'], s=1, color='tab:blue', alpha=0.5)


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

data = FlowCal.io.FCSData('DGY3411.fcs')
channels = data.channels
df = pd.DataFrame(data, columns=channels)
mask = ( df['FSC-A']<4.1e6 ) & ( df['B2-A']>1e4 )
ax.scatter(df[mask]['FSC-A'], df[mask]['B2-A'], s=1, color='tab:olive', alpha=0.5)

data = FlowCal.io.FCSData('DGY3410.fcs')
channels = data.channels
df = pd.DataFrame(data, columns=channels)
mask = ( df['FSC-A']<4.1e6 ) & ( df['B2-A']>1e4 )
ax.scatter(df[mask]['FSC-A'], df[mask]['B2-A'], s=1, color='tab:brown', alpha=0.5)


ax.set_xlabel('Forward scatter area')
ax.set_ylabel('log(mCitrine fluorescence)')
ax.set_yscale('log')
ax.set_ylim(2e4,2e6)
ax.set_title('Copy number from fluorescence')

plt.gcf().align_labels()
plt.tight_layout(h_pad=-0.2)

#plt.savefig('Fig7_combined.png', dpi=500)
#plt.savefig('Fig7_combined.pdf', dpi=500)

plt.show()

#-------------------------------------------------------
