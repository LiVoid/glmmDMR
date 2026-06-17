# Benchmarking Scripts for glmmDMR

This directory contains scripts to compare glmmDMR against other DMR detection methods (DSS, methylKit, etc.) on simulated data with known ground truth.

## Quick Start

Assuming you have already generated simulated data with `../simulate_sites.R`:

```bash
# 1. Convert site-level data to tool-specific formats
Rscript 03.convert_sites_for_otherSoft.R ../results/site_window_sim/tsv/sites_CG.tsv.gz

# 2. Run each comparison tool
Rscript 04.run_DSS.R
Rscript 04.run_methylKit.R

# 3. Run glmmDMR on the same data (outside this directory)
# Rscript ../scripts/run_glmmDMR.R [options...]

# 4. Evaluate results across methods
Rscript 05.evaluate_dmrs.R \
  --simes ../glmmDMR_results/dmrs_simes.tsv \
  --stouffer ../glmmDMR_results/dmrs_stouffer.tsv \
  --combined ../glmmDMR_results/dmrs_combined.tsv \
  --out-prefix results/comparison
```

## Scripts Overview

### 03.convert_sites_for_otherSoft.R

**Purpose:** Convert simulated site-level methylation data into format-specific files for each tool.

**Input:** `sites_CG.tsv.gz` from the simulation

**Output:** 
- `output_for_DSS/sites_CG_forDSS.tsv` (DSS format)
- `output_for_methylKit/sites_CG_forMethylKit_*.txt` (methylKit format, one per sample)

**Usage:**
```bash
Rscript 03.convert_sites_for_otherSoft.R <path/to/sites_CG.tsv.gz>
```

---

### 04.run_DSS.R

**Purpose:** Run DSS (Differentially Methylated Sites via Shrinkage) DMR detection.

**Input:** `output_for_DSS/sites_CG_forDSS.tsv` (created by 03.convert_sites_for_otherSoft.R)

**Output:**
- `output_for_DSS/DSS_dmlTest.tsv` (site-level results)
- `output_for_DSS/DSS_dmrs.tsv` (DMR regions)

**Usage:**
```bash
Rscript 04.run_DSS.R
```

**Parameters (edit in script):**
- `p.threshold`: p-value cutoff for DMR calling (default: 0.05)
- `minlen`: minimum DMR length in bp (default: 500)
- `minCG`: minimum CpG count (default: 3)
- `dis.merge`: merge distance (default: 200 bp)

---

### 04.run_methylKit.R

**Purpose:** Run methylKit DMR detection using tiling + differential methylation.

**Input:** `output_for_methylKit/sites_CG_forMethylKit_*.txt` files (created by 03.convert_sites_for_otherSoft.R)

**Output:**
- `output_for_methylKit/methylKit_diff.tsv` (tile-level differential methylation results)

**Usage:**
```bash
Rscript 04.run_methylKit.R
```

**Parameters (edit in script):**
- `win.size`: tile window size (default: 300 bp)
- `step.size`: tile step size (default: 200 bp)
- `qvalue`: q-value threshold (default: 0.05)
- Coverage filters: `lo.count=5`, `hi.perc=99.9`

---

### 05.evaluate_dmrs.R

**Purpose:** Compare DMR detection results and calculate performance metrics against ground truth.

**Input:** Three DMR result files (from different methods or statistical approaches):
- Simes combined p-value method
- Stouffer combined p-value method
- Combined delta + p-value method

**Output:**
- `<prefix>_distribution.pdf` — DMR size and count distributions
- `<prefix>_performance.pdf` — Sensitivity/specificity/precision curves
- `<prefix>_summary.tsv` — Summary statistics per method

**Usage:**
```bash
Rscript 05.evaluate_dmrs.R \
  --simes <path/to/dmrs_simes.tsv> \
  --stouffer <path/to/dmrs_stouffer.tsv> \
  --combined <path/to/dmrs_combined.tsv> \
  --out-prefix <output/prefix>
```

**Options:**
- `--simes`: Input DMR file from Simes method (required)
- `--stouffer`: Input DMR file from Stouffer method (required)
- `--combined`: Input DMR file from combined method (required)
- `--out-prefix`: Output file prefix (default: `dmr_eval`)

---

## Example: Full Benchmarking Run

```bash
# From simulation/ directory
cd simulation

# Generate simulated data
Rscript simulate_sites.R \
  --out_dir test_sim \
  --seed 12345 \
  --rep_per_group 4 \
  --block_frac 0.05

# Convert formats
cd benchmarking
Rscript 03.convert_sites_for_otherSoft.R ../test_sim/tsv/sites_CG.tsv.gz

# Run DSS
echo "Running DSS..."
time Rscript 04.run_DSS.R

# Run methylKit
echo "Running methylKit..."
time Rscript 04.run_methylKit.R

# (Run glmmDMR in a separate workflow)
# ... glmmDMR results should be in ../glmmDMR_results/

# Evaluate
echo "Evaluating..."
Rscript 05.evaluate_dmrs.R \
  --simes ../glmmDMR_results/windows_CG_fit_beta_pooled_dmr_dmrs_simes.tsv \
  --stouffer ../glmmDMR_results/windows_CG_fit_beta_pooled_dmr_dmrs_stouffer.tsv \
  --combined ../glmmDMR_results/windows_CG_fit_beta_pooled_dmr_dmrs_combined.tsv \
  --out-prefix results/test_comparison

echo "Benchmarking complete!"
ls -lh results/
```

---

## Dependencies

All scripts require R and the following packages:

- **data.table** — Fast data manipulation
- **tidyr** — Data tidying
- **stringr** — String operations
- **ggplot2** — Visualization
- **optparse** — Command-line argument parsing
- **DSS** — For 04.run_DSS.R
- **methylKit** — For 04.run_methylKit.R

### Install missing packages:

```bash
Rscript -e "
  packages <- c('data.table', 'tidyr', 'stringr', 'ggplot2', 'optparse')
  for (pkg in packages) {
    if (!require(pkg, quietly=TRUE)) {
      install.packages(pkg, repos='http://cran.r-project.org')
    }
  }
  if (!require('DSS', quietly=TRUE)) install.packages('DSS')
  if (!require('methylKit', quietly=TRUE)) BiocManager::install('methylKit')
"
```

---

## Notes

- All scripts assume sample groups: `WT01–WT04` (wild-type) and `MT01–MT04` (mutant).
- Modify sample names in scripts if using different group/replicate structures.
- Output directories (e.g., `output_for_DSS/`) are created automatically.
- For reproducibility, document tool versions and random seeds used.

---

## Troubleshooting

**Error: "No such file or directory: output_for_DSS/..."**
- Make sure you ran `03.convert_sites_for_otherSoft.R` first to create the input files.

**Error: "package DSS not found"**
- Install DSS: `Rscript -e "install.packages('DSS')"`

**Error: "sites_CG.tsv.gz not found"**
- Verify that `simulate_sites.R` was run successfully and output is in `../results/site_window_sim/tsv/`

---

## Adding Additional Tools

To benchmark additional DMR detection methods:

1. **Extend 03.convert_sites_for_otherSoft.R** to export tool-specific formats
2. **Create 04.run_TOOLNAME.R** following the pattern of existing scripts
3. **Update 05.evaluate_dmrs.R** to include the new tool in comparisons

See the main [simulation/README.md](../README.md) for more details.
