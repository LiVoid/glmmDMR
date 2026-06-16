<img width="500" height="100" alt="glmmDMR_logo" src=glmmDMR_logo.png />

# glmmDMR

glmmDMR is a statistical framework for replicate-aware detection of differentially methylated regions (DMRs) from whole-genome bisulfite sequencing (WGBS),  enzymatic methyl-seq (EM-seq) or Oxford Nanopore Technologies (ONT) data. The framework integrates generalized linear mixed models (GLMMs) with a seed-based region construction strategy to improve detection accuracy and reduce false positives driven by replicate-level methylation variability.

glmmDMR consists of three core components:

1. **Window-level data organization**: cytosine sites are organized into fixed-length sliding windows while preserving site-level methylation information across samples and replicates.
2. **Replicate-aware statistical modeling**: methylation data within each window are modeled using a GLMM that treats group identity as a fixed effect and biological replicates as random effects, explicitly capturing replicate-level variability without prior aggregation.
3. **Seed-based DMR construction**: window-level statistical signals are integrated into DMRs using one of three seed strategies — `single_seed`, `multi_seed`, or `hybrid_seed` — each initiating region construction from locally high-confidence windows and expanding based on directional consistency and Stouffer-combined p-values.

> **Key finding**: False-positive DMR detections arise predominantly from genomic regions with high replicate-level methylation variability, not from regions with weak biological signals. glmmDMR addresses this directly through GLMM-based variance modeling.

---

## Table of Contents

1. [Repository Contents](#1-repository-contents)
2. [Installation](#2-installation)
3. [Quick Environment Check](#3-quick-environment-check)
4. [Typical Workflow](#4-typical-workflow)
5. [Step-by-Step Inputs and Outputs](#5-step-by-step-inputs-and-outputs)
6. [Runtime Notes](#6-runtime-notes)
7. [Simulation and Benchmarking](#7-simulation-and-benchmarking)
8. [Citation](#8-citation)
9. [License](#9-license)

---

## 1. Repository Contents

**Core scripts:**

| Script | Description |
|---|---|
| `summarize_extractor.py` | Summarize Bismark extractor output (`*.txt.gz`) into per-site methylation counts |
| `modkit_bed_to_binom_input.py` | (Optional, ONT) Convert modkit pileup BED to per-site counts for `BinomTest.py` |
| `BinomTest.py` | Per-site binomial test with bisulfite non-conversion correction and FDR filtering |
| `prepare_matrix.sh` | Build two-group sliding-window matrices from per-site data |
| `run_glmmDMR.R` | Fit GLMM per window; supports `binom`/`beta` × `aggregate`/`site` configurations |
| `merge_window.R` | Integrate significant windows into DMRs using seed-based strategies |
| `make_binned_methylation_bigwig.R` | (Optional) Generate binned methylation bigWig files for visualization |
| `make_binned_variance_bigwig.py` | (Optional) Generate binned variance bigWig files for visualization |

**Simulation resources:**

- `simulation/simulate_sites.R`: simulate synthetic WGBS datasets with known ground-truth DMRs
- `simulation/README.md`: simulation parameter documentation

**Documentation:**

- `tutorial/tutorial_glmmDMR.md`: detailed step-by-step tutorial

---

## 2. Installation

```bash
git clone https://github.com/ktonosaki/glmmDMR.git
cd glmmDMR
```

**Optional: set executable flags**

```bash
chmod +x scripts/summarize_extractor.py scripts/BinomTest.py scripts/prepare_matrix.sh
chmod +x scripts/run_glmmDMR.R scripts/merge_window.R
chmod +x scripts/make_binned_methylation_bigwig.R scripts/make_binned_variance_bigwig.py scripts/modkit_bed_to_binom_input.py
```

**Python dependencies:**

```bash
python -m pip install pandas numpy scipy statsmodels pyBigWig tqdm
```

**R dependencies:**

```r
install.packages(c("optparse", "data.table", "glmmTMB", "future", "future.apply"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "rtracklayer"))
```

**Recommended external tools:**

- [bedtools](https://bedtools.readthedocs.io/)
- [samtools](http://www.htslib.org/)
- gzip / zcat
- [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) (+ Bowtie2) for upstream methylation extraction

---

## 3. Quick Environment Check

```bash
python scripts/summarize_extractor.py --help
python scripts/BinomTest.py --help
python scripts/modkit_bed_to_binom_input.py --help
bash scripts/prepare_matrix.sh --help
Rscript scripts/run_glmmDMR.R --help
Rscript scripts/merge_window.R --help
```

---

## 4. Typical Workflow

```
Raw FASTQ
    ↓  (Bismark alignment + deduplication)
Bismark extractor output
    ↓  summarize_extractor.py      [Step 1]
Per-site methylation counts
    ↓  BinomTest.py                [Step 2]
Per-site binomial test results
    ↓  prepare_matrix.sh           [Step 3]
Sliding-window matrix
    ↓  run_glmmDMR.R               [Step 4]
Window-level GLMM statistics
    ↓  merge_window.R              [Step 5]
DMR calls (TSV + BED)
    ↓  make_binned_*_bigwig        [Step 6, optional]
BigWig tracks for visualization
```

For a detailed tutorial: [tutorial](tutorial/)

---

## 5. Step-by-Step Inputs and Outputs

### Step 1: `summarize_extractor.py`

Summarizes per-context methylation counts from Bismark extractor output files into a single per-site table.

**Input:** Bismark methylation extractor output directory (`*.txt.gz`) for one sample.

**Output:** `*_summarized_output.tsv.gz` with key columns: `chr, pos, strand, meth, unmeth, context`

```bash
python summarize_extractor.py \
  -i sample_methylation_extractor_dir \
  -o sample_summarized_output.tsv.gz \
  --threads 4
```

**Optional (ONT / modkit pileup BED):**

If your input is ONT methylation calls from [modkit](https://github.com/nanoporetech/modkit) `pileup` (BED/BED.GZ), convert it to the same Step-1 output format (`chr, pos, strand, meth, unmeth, context`) before running `BinomTest.py`.

```bash
python scripts/modkit_bed_to_binom_input.py \
  -i sample_modkit_pileup.bed \
  -o sample_summarized_output.tsv.gz \
  --min_coverage 5
```

Then proceed to Step 2 exactly as in the standard workflow.

---

### Step 2: `BinomTest.py`

Performs a per-site binomial test to identify cytosines with methylation levels significantly above the bisulfite non-conversion rate. Non-significant sites are retained in the output but have their `meth` count set to `0`, ensuring they do not contribute to downstream window-level signal while preserving coverage information.

**Input:** `*_summarized_output.tsv.gz`

**Output:** `*_binomtest_result.tsv.gz`

**Estimating the non-conversion rate:**

The null probability (bisulfite non-conversion rate) can be estimated from a non-conversion control chromosome (e.g., chloroplast or lambda spike-in) using `--nonconv_chr`. If no control chromosome is available, provide a precomputed value with `--null_prob`.

```bash
# Recommended: estimate from non-conversion control chromosome
python BinomTest.py \
  -i sample_summarized_output.tsv.gz \
  -o sample_binomtest_result.tsv.gz \
  --nonconv_chr chloroplast \
  --min_coverage 5 \
  --fdr_threshold 0.05 \
  --threads 4
```

```bash
# Alternative: provide precomputed non-conversion rate
python BinomTest.py \
  -i sample_summarized_output.tsv.gz \
  -o sample_binomtest_result.tsv.gz \
  --null_prob 0.012 \
  --min_coverage 5 \
  --fdr_threshold 0.05 \
  --threads 4
```

---

### Step 3: `prepare_matrix.sh`

Builds a two-group sliding-window matrix from per-site binomial test results, organizing methylation counts across all samples and replicates within each window.

**Input:**
- Group 1 and Group 2 sets of `*_binomtest_result.tsv.gz`
- Reference FASTA (with FAI index)

**Output:** Sliding-window matrix TSV.GZ files (one per methylation context)

```bash
bash prepare_matrix.sh \
  --fasta TAIR10.fasta \
  --group1 WT_1_binomtest_result.tsv.gz WT_2_binomtest_result.tsv.gz \
  --group2 MT_1_binomtest_result.tsv.gz MT_2_binomtest_result.tsv.gz \
  --group_labels WT MT \
  --window 500 \
  --slide 300 \
  --output prep_out
```

---

### Step 4: `run_glmmDMR.R`

Fits a GLMM to each window to test for differential methylation between groups while explicitly modeling replicate-level variability as random effects.

**Model configurations** (`--family` × `--mode`):

| Configuration | Response variable | Data representation |
|---|---|---|
| `beta` + `site` | Methylation proportion per cytosine | Site-level (recommended) |
| `beta` + `aggregate` | Window-mean methylation proportion | Aggregated |
| `binom` + `site` | Methylated read counts per cytosine | Site-level |
| `binom` + `aggregate` | Total methylated counts per window | Aggregated |

The **beta site-level** configuration (`--family beta --mode site`) is recommended based on benchmarking: it maintains the highest precision across a wide range of effect sizes and most effectively suppresses false positives driven by high methylation variance.

Group differences are modeled as fixed effects; biological replicates are modeled as random intercepts. Statistical significance is assessed by likelihood ratio test, with p-values adjusted using the Benjamini–Hochberg method.

**Input:** Sliding-window matrix TSV.GZ from `prepare_matrix.sh`

**Output:** `*_fit_<family>_<mode>.tsv.gz` with columns including `chr, start, end, p, delta`

**Key options:**

- `--family`: `binom` or `beta`
- `--mode`: `aggregate` or `site`
- `--random_effect`: enable replicate-level random effects (recommended)
- `--min_cov`: minimum coverage per site (default: 5)
- `--workers`, `--batches`: parallelization settings

```bash
Rscript run_glmmDMR.R \
  -i prep_out/WT_MT_CpG_matrix.tsv.gz \
  -o glmm_out/WT_MT_CpG \
  --group1 WT --group2 MT \
  --family beta \
  --mode site \
  --min_reps_g1 2 --min_reps_g2 2 \
  --min_sites_win 3 --min_cov 5 \
  --random_effect \
  --workers 8 --batches 200
```

---

### Step 5: `merge_window.R`

Integrates window-level GLMM statistics into DMRs using a seed-based region construction strategy. Window-level p-values within candidate regions are combined using the Stouffer method. All windows are first assigned a methylation direction based on the sign of the estimated group effect; only windows sharing the same direction are integrated into a common region.

**Input:** GLMM fit table (`*_fit_<family>_<mode>.tsv.gz`) with required columns: `chr, start, end, p, delta`

**Output:** DMR TSV files and BED files per merge mode

#### Seed strategies

**`single_seed`**: Individual windows with p ≤ `--p-seed` are designated as seeds. Regions are expanded by incorporating adjacent same-direction windows within `--max-gap-bp` bp, provided that inclusion does not increase the Stouffer-combined p-value by more than `--max-p-degradation`-fold.

**`multi_seed`**: Clusters of neighboring same-direction windows are evaluated jointly; a cluster is designated as a seed if its Stouffer-combined p-value falls below `--p-seed`. Expansion follows the same criteria as `single_seed`.

**`hybrid_seed`** (default, recommended): Applies `multi_seed` genome-wide first, then applies `single_seed` to regions not covered by any multi-seed region. This approach balances sensitivity for clustered moderate signals and high-confidence individual signals.

#### Key options

**Seed / Extension:**

| Option | Default | Description |
|---|---|---|
| `--p-seed` | 0.05 | Significance threshold for seed definition |
| `--p-extend` | 0.05 | Significance threshold for region extension |
| `--max-gap-bp` | 200 | Maximum gap (bp) between adjacent windows |
| `--max-p-degradation` | 1.2 | Maximum allowed increase in combined p-value during extension |
| `--max-final-p` | 1.0 | Maximum combined p-value of a retained DMR |
| `--min-strong-windows` | 0.5 | Minimum fraction of windows with p ≤ `--p-seed` |
| `--min-windows` | 1 | Minimum number of windows per DMR |

**Adaptive effect size filter**:

| Option | Default | Description |
|---|---|---|
| `--adaptive-delta` | FALSE | Enable adaptive minimum effect size threshold |
| `--adaptive-delta-method` | `median_ratio` | Method: `median_ratio`, `q50`, `q25`, `q10`, `mad` |
| `--adaptive-delta-ratio` | 0.6 | Ratio for `median_ratio` method |

> **Note:** `--adaptive-delta` performs two passes: an initial pass without an effect size filter to estimate the distribution of Δmethylation across candidate DMRs, followed by a second pass applying a data-driven minimum |Δmethylation| threshold. `q25` is a practical starting point.

**Post-filter:**

| Option | Default | Description |
|---|---|---|
| `--post-filter` | FALSE | Enable post-detection quality filtering |
| `--min-median-p` | 0.01 | Maximum median p-value of windows in a retained DMR |
| `--min-consistent-frac` | 0.5 | Minimum fraction of windows with p ≤ `--p-seed` |

**Overlap merge:**

| Option | Default | Description |
|---|---|---|
| `--merge-overlaps` | FALSE | Re-merge overlapping/nearby same-direction DMRs |
| `--merge-overlaps-gap` | 0 | Maximum gap (bp) for overlap re-merge |

```bash
Rscript merge_window.R \
  --windows glmm_out/WT_MT_CpG_fit_beta_site.tsv.gz \
  --out-prefix dmr_out/WT_MT_CpG \
  --merge-mode hybrid_seed \
  --p-seed 0.05 \
  --p-extend 0.05 \
  --max-p-degradation 1.15 \
  --max-final-p 0.0015 \
  --min-strong-windows 0.5 \
  --min-windows 1 \
  --adaptive-delta
```

---

### Step 6: Optional bigWig tracks

Generate binned methylation level and replicate-level variance bigWig files for genome browser visualization.

**Methylation bigWig** (per sample):

```bash
Rscript make_binned_methylation_bigwig.R \
  -i sample_binomtest_result.tsv.gz \
  -b 50 \
  -o sample_CpG.bw \
  --genome TAIR10.chromInfo \
  --context CpG
```

**Variance bigWig** (across replicates):

```bash
python make_binned_variance_bigwig.py \
  --inputs rep1_CpG.bw rep2_CpG.bw rep3_CpG.bw rep4_CpG.bw \
  --output group_CpG.variance.bw \
  --bin-size 200 \
  --min-tracks 2 \
  --norm none
```

---

## 6. Runtime Notes

- The pipeline assumes two-group comparisons (e.g., wild-type vs. mutant).
- Methylation contexts (CG, CHG, CHH) are processed independently.
- GLMM fitting is the most computationally intensive step. Use `--workers` and `--batches` to tune parallelization and memory usage for your system.
- For large datasets or high-memory jobs, set a dedicated temporary directory:

```bash
export TMPDIR=/path/to/large_storage/tmp
mkdir -p "$TMPDIR"
```

---

## 7. Simulation and Benchmarking

Controlled benchmark datasets with known ground-truth DMRs can be generated using the simulation scripts:

- `simulation/simulate_sites.R`: generates synthetic CG methylation datasets with configurable overdispersion, replicate variability, coverage distribution, and DMR embedding parameters.
- `simulation/README.md`: full documentation of simulation parameters.

For details on benchmarking methodology and comparison with existing DMR detection methods (DSS, methylKit, DMRfinder, metilene, MACAU2, Fisher's exact test), see the associated publication.

---

## 8. Citation

If you use glmmDMR in your research, please cite:

> [Author names] (2025). glmmDMR: a GLMM-based framework for replicate-aware detection of differentially methylated regions. ****. [DOI upon publication]

---

## 9. License

glmmDMR is released under the MIT License. See [LICENSE](LICENSE) for details.

This software relies on third-party R packages including glmmTMB, GenomicRanges, and others listed in the Installation section. Please refer to the respective package documentation for their licensing terms.
