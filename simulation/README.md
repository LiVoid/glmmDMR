# Create simulate data set (CG context)

This tutorial explains what `simulate_sites.R` does, how to run it, what each option means, and what files are produced.

The script simulates site-level methylation count data for two groups (for example, `WT` vs `MT`) with known DMR truth labels.

## 1. What This Script Generates

For a fixed context (`CG`), the script simulates:

1. Genomic DMR blocks with random directions (`hyper` or `hypo`) and effect sizes (`delta`).
2. Site positions across a synthetic chromosome.
3. Group- and replicate-specific methylation counts (`meth`, `unmeth`) with:
   - beta-binomial overdispersion,
   - replicate effects,
   - optional site-level heterogeneity,
   - missingness and coverage filtering.
4. Ground-truth labels for both blocks and sites.

## 2. Quick Start

Run with defaults:

```bash
Rscript simulate_sites.R
```

Run with custom settings:

```bash
Rscript simulate_sites.R \
  --out_dir results/site_window_sim \
  --chr_len 10000000 \
  --site_lambda 20000 \
  --block_frac 0.05 \
  --block_len_lo 300 \
  --block_len_hi 2000 \
  --delta_grid 0.1,0.2,0.3,0.4 \
  --groups WT,MT \
  --rep_per_group 4 \
  --base 0.7 \
  --dmr_site_prop 0.5 \
  --mu_site_sd 0.05 \
  --rho 0.3 \
  --logcov_mu 3.0 \
  --logcov_sd 1.2 \
  --miss_rate 0.30 \
  --min_cov 5 \
  --seed 202509
```

## 3. Model Summary

The script uses a CG-only simulation model.

- Base methylation level is set by `--base`.
- DMR effect (`delta`) is applied to group 2 only.
- If a site overlaps multiple DMR blocks, the delta with the largest absolute value is used.
- Per-site probability is sampled via beta distribution (with concentration derived from `rho`), then counts are drawn from binomial.

Informally:

```text
mu = base + delta * I(group == second_group) + replicate_effect + optional_site_effect
p ~ Beta(mu, concentration)
meth ~ Binomial(coverage, p)
unmeth = coverage - meth
```

## 4. Option Reference

## Required/General

`--out_dir`
- Output directory root.
- Default: `results/site_window_sim`

`--seed`
- Random seed for reproducibility.
- Default: `202509`

## Genome/Site Layout

`--chr_len`
- Simulated chromosome length in bp.
- Default: `5e6`

`--site_lambda`
- Site density per Mb for CG.
- Default: `4e4`

## DMR Block Construction

`--block_frac`
- Approximate fraction of chromosome covered by candidate DMR blocks.
- Default: `0.10`

`--block_len_lo`
- Minimum DMR block length (bp).
- Default: `500`

`--block_len_hi`
- Maximum DMR block length (bp).
- Default: `5000`

`--delta_grid`
- Comma-separated absolute effect sizes to sample from.
- Sign is assigned randomly (`hyper`/`hypo`).
- Default: `0.1,0.2,0.3`

## Group/Replicate Design

`--groups`
- Two group labels, comma-separated.
- Must contain exactly 2 groups.
- Default: `WT,MT`

`--rep_per_group`
- Number of replicates per group.
- Default: `3`

## Methylation Mean Structure

`--base`
- Baseline methylation rate for CG in group 1.
- Default: `0.7`

`--dmr_site_prop`
- Fraction of overlapping sites that receive DMR effect inside blocks.
- Range: `[0, 1]`
- Default: `1.0`

`--mu_site_sd`
- SD of optional site-level random effect added to mean methylation.
- Default: `0.0`

## Overdispersion and Coverage

`--rho`
- Base overdispersion scale for beta-binomial component.
- Default: `0.05`

`--logcov_mu`
- Mean of log-coverage distribution.
- Default: `log(20)`

`--logcov_sd`
- SD of log-coverage distribution.
- Default: `0.5`

`--miss_rate`
- Probability that a site is set to missing (`coverage = 0`) before coverage filtering.
- Default: `0.10`

`--min_cov`
- Minimum coverage threshold to keep a site.
- Default: `5`

## 5. Outputs

All outputs are written to:

```text
<out_dir>/tsv/
```

with file names:

1. `sites_CG.tsv.gz`
2. `truth_blocks_CG.tsv.gz`
3. `truth_sites_CG.tsv.gz`

### 5.1 sites_CG.tsv.gz

Site-level simulated observations with truth labels joined.

Columns:
- `chr`: chromosome name (`chr1`)
- `pos`: genomic position (1-based)
- `sample`: sample ID (for example `WT01`)
- `group`: group label
- `replicate`: replicate index
- `context`: fixed as `CG`
- `meth`: methylated read count
- `unmeth`: unmethylated read count
- `truth`: `1` if site overlaps any truth block, else `0`
- `dir`: `hyper`, `hypo`, or `NA`

### 5.2 truth_blocks_CG.tsv.gz

Truth DMR blocks used during simulation.

Columns:
- `chr`
- `start`
- `end`
- `dir`
- `delta`

### 5.3 truth_sites_CG.tsv.gz

Unique sites with truth labels.

Columns:
- `chr`
- `pos`
- `truth`
- `dir`

## 6. Practical Notes

- The script is currently CG-only (`ctx <- "CG"`).
- `--groups` must provide exactly two labels.
- `--dmr_site_prop` must be in `[0, 1]`.
- For overlapping DMR blocks at one site, the script keeps the effect with largest absolute `delta`.

## 7. Suggested Parameter Sets

More realistic / moderate noise:

```bash
--site_lambda 20000 --block_frac 0.05 --rho 0.10 --miss_rate 0.10 --dmr_site_prop 0.8
```

Hard setting (more missingness/noise):

```bash
--site_lambda 20000 --block_frac 0.05 --rho 0.30 --miss_rate 0.30 --dmr_site_prop 0.5 --mu_site_sd 0.05
```

## 8. Minimal Validation After Run

```bash
zcat results/site_window_sim/tsv/sites_CG.tsv.gz | head
zcat results/site_window_sim/tsv/truth_blocks_CG.tsv.gz | head
zcat results/site_window_sim/tsv/truth_sites_CG.tsv.gz | head
```

---

# Benchmarking: Compare glmmDMR Against Other Tools

This section describes how to benchmark glmmDMR against other commonly used DMR detection methods (DSS, methylKit, Fisher's exact test, etc.) using simulated data with known ground truth.

## Overview

The benchmarking workflow has three main stages:

1. **Data Format Conversion** → Convert simulated site-level data into tool-specific formats
2. **Run Comparison Tools** → Execute each DMR detection method independently
3. **Evaluate Results** → Compare performance (sensitivity, specificity, runtime, etc.)

All scripts are located in the `benchmarking/` subdirectory.

## Prerequisites

Ensure that **all required R packages for each tool are installed**:

```bash
# DSS and methylKit (if not already installed)
Rscript -e "
  if (!require('DSS', quietly=TRUE)) install.packages('DSS', repos='http://cran.r-project.org')
  if (!require('methylKit', quietly=TRUE)) BiocManager::install('methylKit')
"
```

## Stage 1: Format Conversion

Convert the simulated site-level data into the input formats required by each comparison tool.

### Script: `benchmarking/03.convert_sites_for_otherSoft.R`

**Purpose:** Convert `sites_CG.tsv.gz` into tool-specific formats.

**Usage:**

```bash
cd benchmarking
Rscript 03.convert_sites_for_otherSoft.R ../results/site_window_sim/tsv/sites_CG.tsv.gz
```

**Inputs:**
- `sites_CG.tsv.gz` from simulation (columns: chr, pos, sample, group, replicate, context, meth, unmeth, truth, dir)

**Outputs:**

The script creates two directories:

- `output_for_DSS/`
  - `sites_CG_forDSS.tsv` — Format for DSS (columns: chr, pos, N, X, sample)

- `output_for_methylKit/`
  - `sites_CG_forMethylKit_{MT1,MT2,MT3,MT4,WT1,WT2,WT3,WT4}.txt` — Per-sample methylKit format
  - Columns: chrBase, chr, base, strand, coverage, freqC, freqT

**Note:** The script includes conversion logic for metilene and other tools as well; extend as needed.

---

## Stage 2: Run Comparison Tools

Each tool is run independently on the converted data.

### DSS (Differentially Methylated Sites via Shrinkage)

**Script:** `benchmarking/04.run_DSS.R`

**Usage:**

```bash
cd benchmarking
Rscript 04.run_DSS.R
```

**What it does:**
- Reads DSS-formatted data from `output_for_DSS/sites_CG_forDSS.tsv`
- Runs `DMLtest()` and `callDMR()` with default parameters:
  - Smoothing span: 500 bp
  - P-value threshold: 0.05
  - Minimum DMR length: 500 bp
  - Minimum CpG count: 3

**Outputs:**
- `output_for_DSS/DSS_dmlTest.tsv` — Site-level test results
- `output_for_DSS/DSS_dmrs.tsv` — Detected DMRs

**Customization:** Edit script to adjust `p.threshold`, `minlen`, `minCG`, `dis.merge` parameters.

---

### methylKit (Diffmeth + Tiling)

**Script:** `benchmarking/04.run_methylKit.R`

**Usage:**

```bash
cd benchmarking
Rscript 04.run_methylKit.R
```

**What it does:**
- Reads per-sample methylKit files from `output_for_methylKit/`
- Filters coverage: `lo.count=5`, `hi.perc=99.9`
- Creates tiles (window size 300 bp, step 200 bp)
- Calls differential methylation with `qvalue=0.05`

**Outputs:**
- `output_for_methylKit/methylKit_diff.tsv` — Tile-level results

**Customization:** Edit script to adjust tile size (`win.size`), step (`step.size`), or q-value threshold.

---

## Stage 3: Evaluate Results

Compare detected DMRs across all methods against ground truth.

### Script: `benchmarking/05.evaluate_dmrs.R`

**Purpose:** Compare DMR detection methods (e.g., glmmDMR Simes vs Stouffer vs combined) and visualize performance metrics.

**Usage:**

```bash
cd benchmarking
Rscript 05.evaluate_dmrs.R \
  --simes ../glmmDMR_results/dmrs_simes.tsv \
  --stouffer ../glmmDMR_results/dmrs_stouffer.tsv \
  --combined ../glmmDMR_results/dmrs_combined.tsv \
  --out-prefix ../evaluation_results/dmr_eval
```

**Inputs:**
- Three DMR result files (TSV format) with columns: chr, start, end, [other columns]
- Each file represents one DMR detection method or statistical approach

**Outputs:**
- `dmr_eval_distribution.pdf` — DMR length/count distributions by method
- `dmr_eval_performance.pdf` — Sensitivity/specificity/precision curves
- `dmr_eval_summary.tsv` — Summary statistics (count, median length, etc.)

**Comparison Metrics:**
- **Sensitivity:** Fraction of true DMRs detected
- **Specificity:** Fraction of non-DMRs correctly excluded
- **Precision:** Fraction of detected regions that overlap truth

---

## Complete Benchmarking Workflow Example

```bash
# Assuming you have run: Rscript simulate_sites.R

# Step 1: Convert formats
cd benchmarking
Rscript 03.convert_sites_for_otherSoft.R ../results/site_window_sim/tsv/sites_CG.tsv.gz

# Step 2: Run each tool
echo "Running DSS..."
time Rscript 04.run_DSS.R

echo "Running methylKit..."
time Rscript 04.run_methylKit.R

# (Also run glmmDMR on the same data using ../scripts/run_glmmDMR.R)

# Step 3: Evaluate results
echo "Evaluating results..."
Rscript 05.evaluate_dmrs.R \
  --simes ../glmmDMR_results/windows_CG_fit_beta_pooled_dmr_dmrs_simes.tsv \
  --stouffer ../glmmDMR_results/windows_CG_fit_beta_pooled_dmr_dmrs_stouffer.tsv \
  --combined ../glmmDMR_results/windows_CG_fit_beta_pooled_dmr_dmrs_combined.tsv \
  --out-prefix results/comparison

echo "Done! Results in results/comparison_*.pdf and results/comparison_summary.tsv"
```

---

## Adding New Tools

To integrate a new DMR detection tool:

1. **Add format conversion** in `03.convert_sites_for_otherSoft.R` (create new `output_for_TOOL/` directory)
2. **Create run script** `04.run_TOOL.R` with:
   - Tool-specific parameter setup
   - DMR calling code
   - Output saved as TSV (columns: chr, start, end, [p-value, etc.])
3. **Update evaluation script** to include the new tool's results for comparison

---

## Notes

- All scripts assume the same sample labels and group structure: `WT01–WT04` (group 1) and `MT01–MT04` (group 2).
- Paths to input/output directories can be customized by editing the scripts.
- For reproducible benchmarking, ensure all random seeds are fixed and tool versions are documented.

