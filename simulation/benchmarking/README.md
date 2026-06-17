# Benchmarking scripts for glmmDMR

This directory contains benchmarking scripts used to compare glmmDMR with DSS, methylKit, Fisher, metilene, DMRfinder, and MACAU2.

This README mirrors the practical workflow used in bash_variance.sh and clarifies which extra processing steps are done inside wrapper scripts.

## Files in this directory

- 03.convert_sites_for_otherSoft.R
- run_DSS.R
- run_methylKit.R
- run_fisher.R
- 04.run_metilene.R
- 04.run_dmrfinder.R
- 04.run_MACAU2.R
- findDMRs_fixed.r
- evaluate_dmrs.R

## Important compatibility notes

- Use the file names above exactly. Some older notes refer to 04.run_DSS.R, 04.run_methylKit.R, and 05.evaluate_dmrs.R, but this directory uses run_DSS.R, run_methylKit.R, and evaluate_dmrs.R.
- 03.convert_sites_for_otherSoft.R requires one input argument: the sites file. If extra positional arguments are provided, they are ignored.
- MACAU2 is site-level analysis. The script also writes merged significant windows as a convenience output.

## External tools and scripts required

These are not bundled in this repository:

- metilene binary
- metilene_output.pl
- DMRfinder combine_CpG_sites.py
- MACAU2 package or MACAU2 R source directory

## R package dependencies

- data.table
- tidyr
- stringr
- ggplot2
- optparse
- dplyr
- readr
- future
- furrr
- DSS
- methylKit
- MACAU2

## Install R dependencies

```bash
Rscript -e "
packages <- c('data.table','tidyr','stringr','ggplot2','optparse','dplyr','readr','future','furrr')
for (pkg in packages) {
  if (!require(pkg, quietly=TRUE)) install.packages(pkg, repos='https://cloud.r-project.org')
}
if (!require('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org')
if (!require('DSS', quietly=TRUE)) BiocManager::install('DSS')
if (!require('methylKit', quietly=TRUE)) BiocManager::install('methylKit')
"
```

## Workflow

### 1. Convert site data to tool-specific inputs

```bash
Rscript 03.convert_sites_for_otherSoft.R ../results/site_window_sim/tsv/sites_CG.tsv.gz
```

This creates:

- output_for_DSS/sites_CG_forDSS.tsv
- output_for_methylKit/sites_CG_forMethylKit_*.txt
- output_for_metilene/sites_CG_forMetilene.txt
- output_for_DMRfinder/sites_CG_forDMRfinder_*.txt
- output_for_MACAU/windows_CG_forMACAU_*.bed

### 2. Run each method

DSS:

```bash
Rscript run_DSS.R
```

methylKit:

```bash
Rscript run_methylKit.R
```

Fisher (window input):

```bash
Rscript run_fisher.R \
  -i ../results/site_window_sim/tsv/windows_CG.tsv.gz \
  --input_type windows \
  -c CG \
  -o output_for_fisher/fisher_out.tsv \
  --merge_dmrs \
  --fdr_threshold 0.05 \
  --max_gap_bp 300 \
  --min_windows 1 \
  --threads 36
```

metilene:

```bash
Rscript 04.run_metilene.R \
  --metilene-bin /path/to/metilene \
  --metilene-output-pl /path/to/metilene_output.pl
```

DMRfinder:

```bash
Rscript 04.run_dmrfinder.R \
  --python-bin python \
  --combine-script /path/to/combine_CpG_sites.py \
  --finddmrs-script ./findDMRs_fixed.r
```

MACAU2:

```bash
Rscript 04.run_MACAU2.R \
  --macau2-r-dir /path/to/MACAU2/R
```

### 3. Evaluate glmmDMR outputs

```bash
Rscript evaluate_dmrs.R \
  --simes ../glmmDMR_results/dmrs_simes.tsv \
  --stouffer ../glmmDMR_results/dmrs_stouffer.tsv \
  --combined ../glmmDMR_results/dmrs_combined.tsv \
  --out-prefix results/comparison
```

## What is additionally processed in wrapper scripts

Compared with direct bash lines, wrappers handle these processing details:

- 04.run_dmrfinder.R
  - Runs combine_CpG_sites.py
  - Applies sample-name cleanup in results.mod.csv
  - Appends one dummy row for chr1-only edge cases
  - Runs findDMRs_fixed.r

- 04.run_metilene.R
  - Runs metilene and writes metilene_out.tsv
  - Runs metilene_output.pl and writes filter output

- 04.run_MACAU2.R
  - Builds count and coverage matrices from per-sample BED files
  - Runs MACAU2 in BMM mode
  - Writes site-level results and merged significant windows

## Quick output check

```bash
head output_for_DSS/DSS_dmrs.tsv
head output_for_methylKit/methylKit_diff.tsv
head output_for_metilene/metilene_out.tsv
head output_for_DMRfinder/out_findDMRs.txt
head output_for_fisher/fisher_out.tsv
head output_for_MACAU/MACAU2_sites.tsv
```

## Troubleshooting

- Missing output_for_* files
  - Run 03.convert_sites_for_otherSoft.R first.

- metilene fails
  - Check metilene path and metilene_output.pl path.
  - Check that the input matrix has expected sample columns.

- DMRfinder fails
  - Verify combine_CpG_sites.py path.
  - Verify findDMRs_fixed.r exists in this directory or pass --finddmrs-script.
  - Verify Python and Rscript are available.

- MACAU2 fails
  - Install MACAU2 package or pass --macau2-r-dir.
  - Check that all windows_CG_forMACAU_*.bed files exist.
