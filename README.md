# CNV Stability Without Selection Pressure

Code and data associated with **De et al. 2025**. Each figure script contains all required data loading, analysis, and visualization for the corresponding figure in the paper.  

**Preprint:** https://doi.org/10.1101/2025.07.21.665951

Claude AI was used only to write the README below. All code development, data generation and analysis were done by Titir De. 

---

## Repository Structure

```
.
├── Fig1.py                          # Figure 1: CNV strain schematics
├── Fig2.py                          # Figure 2: Relative fitness of CNV strains
├── Fig3.py                          # Figure 3: CNV dynamics across LTEE lineages
├── Fig4.py                          # Figure 4: Phases of CNV loss and time-to-loss
├── Fig6.py                          # Figure 6: Fitness of revertants vs. CNV parents
├── Fig7.py                          # Figure 7: PacBio read depth and flow cytometry
└── Data/
    ├── Fig_Data/                    # Input data files for all figures
    │   ├── *.csv                    # Tabular data (see below)
    │   ├── *.fcs                    # Flow cytometry files (Figure 7)
    │   └── Bedgraph_files           # Link to externally hosted bedgraph/SAM files
    └── Read_Depth_Plots_Revertants/ # Read-depth plots for revertant strains (PNG)
```

---

## Figure Scripts

### `Fig1.py` — CNV Strain Schematics

Produces two panels:

- **Fig 1A:** Schematic illustration of the three possible outcomes for CNV frequency over time — selection for CNVs, stable CNVs, or loss of CNVs — drawn as sigmoid/flat/sigmoidal curves with added noise.
- **Fig 1B:** Genomic maps of all CNV strains used in the study, drawn as colored filled regions on chromosomes XI (GAP1 locus), XIV (MEP2 locus), and XV (PUT4 locus). A fourth panel shows double-aneuploid strains carrying CNVs on both Chr XI and Chr XIV.

**Key inputs:** `All_strain_schematic_GAP1_20250122.CSV`, `All_strain_schematic_MEP2.CSV`, `All_strain_schematic_PUT4.CSV`

---

### `Fig2.py` — Relative Fitness of CNV Strains

Produces three panels:

- **Fig 2A:** Raw CNV frequency (percentage of cells with CNV) across 15 generations of pairwise competition for three replicates of strain GM_13 (DGY2750).
- **Fig 2B:** Log-ratio of strain to ancestor frequency over time, with a linear regression fit and 95% confidence interval. The slope estimates relative fitness.
- **Fig 2C:** Scatter plot of CNV fitness cost (calculated from gene expression data) vs. measured relative fitness across all CNV strains, with a Monte Carlo regression and 95% CI band. Tests whether larger/costlier CNVs are less fit.
- **Fig 2D:** Bar chart of relative fitness (slope + 1) for all CNV strains, with error bars (2 × SE) and significance stars for strains with p ≤ 0.05.

**Statistical methods:** Linear regression (`scipy.stats.linregress`), Monte Carlo simulation (10,000 iterations) for regression uncertainty, Pearson correlation.

**Key inputs:** `Pairwise_III_edited_new.CSV`, `GaschFitnessCost_vs_Size.CSV`, `Pairwise_Regression_Norm_All_30C_final.CSV`

---

### `Fig3.py` — CNV Dynamics Across LTEE Lineages

Produces three multi-panel figures tracking the percentage of cells carrying a CNV over hundreds of generations in a long-term evolution experiment (LTEE):

- **Fig 3A:** Four-panel plot showing CNV dynamics for three GAP1-CNV strains (G_1, G_2, G_3) and control strains (multi-copy and one-copy). Each main panel shows LTEE data; each inset shows the ancestral evolution experiment from which the strain was isolated.
- **Fig 3B:** Eight-panel plot (2 × 4 grid) showing dynamics for eight additional strains: four GAP1-CNV strains (G_4–G_9), one MEP2-CNV strain (M_10), and one PUT4-CNV strain (P_11).
- **Fig 3C:** Eight-panel plot (2 × 4 grid) for additional MEP2 and PUT4 CNV strains plus double aneuploids (strains GM_12–GM_15).

**Key inputs:** `2024_11_Percentage_gating_LTEE3_YPD_20241205_edited.csv`, `2023_12_Percentage_gating_LTEE2_edit.csv`, `2018_Percentage_CNV_acquire_Lauer.csv`, `2018_Percentage_CNV_acquire_Rahman.csv`, `2018_Percentage_CNV_acquire_Avecilla.csv`

---

### `Fig4.py` — Phases of CNV Loss and Time-to-Loss

Produces a four-panel figure examining when CNV loss occurs during the LTEE:

- **Fig 4A:** Example trajectory (DGY1886_1) with a cubic spline fit, annotated with the 75%, 50%, and 25% CNV-frequency thresholds that define the early, middle, and late phases of CNV loss.
- **Figs 4B–4D:** Box plots comparing the generation at which CNV frequency drops below 75% (early), 50% (middle), and 25% (late) for three CNV classes — Chr XI segmental (G_5), Chr XI aneuploid (GM strains), and Chr XIV aneuploid (GM strains). Points exceeding 220 generations are plotted as triangles above the y-axis; sample sizes are annotated. Significance brackets show Welch's t-test p-values between groups.

**Statistical methods:** Cubic spline interpolation (`scipy.interpolate.CubicSpline`), root-finding (`scipy.optimize.root_scalar`), Welch's t-test (`scipy.stats.ttest_ind`), pairwise comparisons across all three CNV classes.

**Key inputs:** `2023_12_Percentage_gating_LTEE2_edit.csv`

---

### `Fig6.py` — Fitness of Revertants vs. CNV Parents

Produces a single bar chart comparing relative fitness (slope + 1 from pairwise competition assays) for all CNV parent strains and their revertants (both complete and partial):

- Bars are colored by type: grey for CNV strains, light grey for complete revertants, hatched for partial revertants.
- Error bars show 2× SE; dotted vertical lines separate strain groups (G_4, GM_12–GM_15).
- Cyan shaded rectangles on CNV-parent bars show fitness ranges estimated by simulation-based inference (SBI).

**Key inputs:** `Pairwise_Regression_Norm_All_30C_final.CSV`

---

### `Fig7.py` — PacBio Read Depth and Flow Cytometry

Produces a figure combining genomic and phenotypic data from a CNV parent strain (DGY1886) and its revertants (DGY3417, DGY3410, DGY3414, DGY3418, DGY3420):

- **Left panels (ax1–ax3):** PacBio read depth across Chr XI (positions 400–650 kb), normalized to ploidy, with a Gaussian smoothing filter (σ = 20). Shows the CNV parent (ax1), partial revertant (ax2), and complete revertants (ax3, four strains overlaid). Vertical dotted lines mark CNV breakpoints; split reads (from SAM files) indicating structural variant junctions are plotted as red tick marks on ax1 and ax2.
- **Right panel (ax4):** Flow cytometry scatter plot (FSC-A vs. log mCitrine fluorescence) for all six strains, used to confirm CNV copy number by fluorescence intensity. Data are loaded from `.fcs` files using `FlowCal`.

**External data:** Bedgraph and SAM files are large and hosted externally; see `Data/Fig_Data/Bedgraph_files` for the download link.

**Key inputs:** `PacBio_*.bedgraph` and `Split_reads_*.sam` (external), `DGY*.fcs`

---

## Data Files

All tabular data files are in `Data/Fig_Data/`:

| File | Description | Used in |
|------|-------------|---------|
| `All_strain_schematic_GAP1_20250122.CSV` | Genomic coordinates (start/end) and copy number for all GAP1-CNV strains on Chr XI | Fig 1 |
| `All_strain_schematic_MEP2.CSV` | Genomic coordinates for MEP2-CNV strains on Chr XIV | Fig 1 |
| `All_strain_schematic_PUT4.CSV` | Genomic coordinates for PUT4-CNV strains on Chr XV | Fig 1 |
| `Pairwise_III_edited_new.CSV` | Pairwise competition data (% CNV by generation) for strain DGY2750 replicates | Fig 2 |
| `GaschFitnessCost_vs_Size.CSV` | Calculated CNV fitness costs and measured fitness slopes/SE for all strains | Fig 2 |
| `Pairwise_Regression_Norm_All_30C_final.CSV` | Linear regression results (slope, intercept, SE, p-value, type) for all strains and revertants | Figs 2, 6 |
| `2024_11_Percentage_gating_LTEE3_YPD_20241205_edited.csv` | Flow cytometry gating data (% multicopy mCitrine/mCherry) for LTEE3 strains by generation | Fig 3 |
| `2023_12_Percentage_gating_LTEE2_edit.csv` | Flow cytometry gating data for LTEE2 strains by generation | Figs 3, 4 |
| `2018_Percentage_CNV_acquire_Lauer.csv` | Historical CNV acquisition data (Lauer et al. 2018) used for inset plots | Fig 3 |
| `2018_Percentage_CNV_acquire_Rahman.csv` | Historical CNV acquisition data (Rahman et al. 2025) used for inset plots | Fig 3 |
| `2018_Percentage_CNV_acquire_Avecilla.csv` | Historical CNV acquisition data (Avecilla et al. 2022) used for inset plots | Fig 3 |
| `Bedgraph_files` | Text file containing a Google Drive link to PacBio bedgraph and split-read SAM files | Fig 7 |
| `DGY*.fcs` | Flow cytometry FCS files for CNV parent and revertant strains | Fig 7 |

`Data/Read_Depth_Plots_Revertants/` contains read-depth PNG plots for individual revertant strains on Chr XI and Chr XIV (named `Chr11_<DGY#>.png` and `Chr14_<DGY#>.png`).

---

## Dependencies

```
numpy
matplotlib
pandas
scipy
statistics (standard library)
FlowCal        # Fig7 only — flow cytometry FCS file parsing
```

Install with:
```bash
pip install numpy matplotlib pandas scipy FlowCal
```

`pysam` is imported in `Fig7.py` but not required (SAM files are parsed directly with `pandas`).

---

## Usage

Each script is self-contained. Run from the directory containing the data files, or update the file paths at the top of each script. Output is displayed interactively; commented-out `plt.savefig()` calls at the bottom of each script can be uncommented to save figures to disk.
