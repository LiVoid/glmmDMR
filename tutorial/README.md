# glmmDMR Tutorial (Detailed Reference)

This document is a practical, detailed reference for glmmDMR, covering what each script does, key options, and input/output formats. For a quick overview of the pipeline, see [README.md](../README.md).

**Target scripts:**
- `summarize_extractor.py`
- `BinomTest.py`
- `prepare_matrix.sh`
- `run_glmmDMR.R`
- `merge_window.R`
- `make_binned_methylation_bigwig.R` (optional)
- `make_binned_variance_bigwig.py` (optional)

---

## 1. Pipeline Summary

**Standard workflow:**

```
Raw FASTQ
    ↓  Quality trimming (fastp)                [Upstream, Section 2]
    ↓  Bismark alignment + deduplication       [Upstream, Section 2]
Bismark extractor output
    ↓  summarize_extractor.py                  [Step 1, Section 3.1]
Per-site methylation counts
    ↓  BinomTest.py                            [Step 2, Section 3.2]  ← per sample
Per-site binomial test results
    ↓  prepare_matrix.sh                       [Step 3, Section 3.3]  ← integrates all samples
Sliding-window matrix
    ↓  run_glmmDMR.R                           [Step 4, Section 3.4]
Window-level GLMM statistics
    ↓  merge_window.R                          [Step 5, Section 3.5]
DMR calls (TSV + BED)
    ↓  make_binned_*_bigwig                    [Step 6, Section 3.6-3.7, optional]
BigWig tracks for visualization
```

**Execution unit guide:**
- Steps 1–2 (Sections 3.1–3.2) are run **per sample**.
- Step 3 (Section 3.3) integrates all per-sample outputs by group.
- Steps 4–5 operate on the integrated window matrix.

---

## 2. Upstream Analysis (before glmmDMR)

Below is a minimal example of the upstream pipeline used to generate Bismark extractor output from raw FASTQ files.

**Quality trimming:**

```bash
fastp \
  -i sample_R1.fastq.gz \
  -I sample_R2.fastq.gz \
  -o clean/sample_R1.trimmed.fq.gz \
  -O clean/sample_R2.trimmed.fq.gz \
  --qualified_quality_phred 30 \
  --cut_front --cut_tail \
  --cut_window_size 4 \
  --cut_mean_quality 20 \
  --length_required 36 \
  --thread 8
```

**Alignment:**

```bash
bismark \
  --bowtie2 \
  -p 8 \
  -o align \
  /path/to/bismark_genome \
  -1 clean/sample_R1.trimmed.fq.gz \
  -2 clean/sample_R2.trimmed.fq.gz

samtools view -@ 8 -q 42 -b align/sample_R1.trimmed_bismark_bt2_pe.bam | \
samtools sort -n -@ 8 -o align/sample.Q42.bam
```

> **Note:** Alignments with mapping quality (MAPQ) below 42 are discarded to reduce multi-mapping artifacts.

**Deduplication and methylation calling:**

```bash
# PCR deduplication
deduplicate_bismark \
  --paired \
  --output_dir align \
  --bam \
  align/sample.Q42.bam

# Methylation calling
bismark_methylation_extractor \
  --paired-end \
  --output calls/sample \
  --parallel 8 \
  --buffer_size 30% \
  --gzip \
  --genome_folder /path/to/bismark_genome \
  align/sample.Q42.deduplicated.bam
```

Use the resulting output files (`CpG_*.txt.gz`, `CHG_*.txt.gz`, `CHH_*.txt.gz`) as input to `summarize_extractor.py` in Step 1.

---

## 3. Script-by-Script Details

### 3.1 `summarize_extractor.py`

Scans Bismark extractor output files (`CpG_*.txt.gz`, `CHG_*.txt.gz`, `CHH_*.txt.gz`) in the input directory, aggregates methylated and unmethylated counts per cytosine site, and writes a merged per-site table across all three contexts.

**Options:**

| Option | Default | Description |
|---|---|---|
| `-i, --input` | required | Bismark extractor output directory |
| `-o, --output` | required | Output `*.tsv.gz` |
| `--threads` | 4 | Number of worker threads |
| `--keep-temp` | FALSE | Keep context-specific temporary files |

**Input format:** `*.txt.gz` files containing Bismark extractor 5-column records (read_id, strand, chr, pos, code).

**Output format (`*_summarized_output.tsv.gz`):**

| Column | Type | Description |
|---|---|---|
| `chr` | string | Chromosome name |
| `pos` | int (1-based) | Genomic position |
| `strand` | `+` / `-` | Strand |
| `meth` | int | Methylated read count |
| `unmeth` | int | Unmethylated read count |
| `context` | `CpG` / `CHG` / `CHH` | Cytosine context |

**Example:**

```bash
python summarize_extractor.py \
  -i calls/sample \
  -o matrix/sample_summarized_output.tsv.gz \
  --threads 4
```

---

### 3.2 `BinomTest.py`

Performs a per-site binomial test to identify cytosines with methylation levels significantly above the bisulfite non-conversion rate. The test is applied independently per context, with p-values adjusted using the Benjamini–Hochberg method.

> **Behavior note:** Non-significant sites are **retained** in the output but have their `meth` count set to `0`. This ensures that such sites do not contribute to downstream window-level signal while preserving coverage information for accurate replicate filtering and window construction in `prepare_matrix.sh`.

**Options:**

| Option | Default | Description |
|---|---|---|
| `-i, --input` | required | `summarize_extractor.py` output |
| `-o, --output` | required | Output `*.tsv.gz` |
| `--nonconv_chr` | None | Chromosome used to estimate non-conversion rate (recommended) |
| `--null_prob` | None | Precomputed non-conversion rate (use when `--nonconv_chr` is unavailable) |
| `--fdr_threshold` | 0.05 | FDR threshold |
| `--min_coverage` | 0 | Minimum coverage per site |
| `--threads` | 4 | Number of worker threads |

**Estimating the non-conversion rate:**

- **Recommended:** Use `--nonconv_chr` to estimate the null probability from a non-conversion control chromosome (e.g., chloroplast, lambda spike-in). This is the most accurate approach.
- **Alternative:** If no control chromosome is available, provide a precomputed value via `--null_prob`. A typical value for WGBS of plants is 0.001–0.015, depending on the bisulfite conversion efficiency.

**Output format (`*_binomtest_result.tsv.gz`):**

| Column | Type | Description |
|---|---|---|
| `chr` | string | Chromosome name |
| `pos` | int (1-based) | Genomic position |
| `strand` | `+` / `-` | Strand |
| `meth` | int | Methylated count (set to 0 for non-significant sites) |
| `unmeth` | int | Unmethylated count |
| `context` | `CpG` / `CHG` / `CHH` | Cytosine context |

**Example (with non-conversion control chromosome):**

```bash
python BinomTest.py \
  -i matrix/sample_summarized_output.tsv.gz \
  -o binom/sample_binomtest_result.tsv.gz \
  --nonconv_chr chloroplast \
  --min_coverage 5 \
  --fdr_threshold 0.05 \
  --threads 4
```

**Example (without non-conversion control chromosome):**

```bash
python BinomTest.py \
  -i matrix/sample_summarized_output.tsv.gz \
  -o binom/sample_binomtest_result.tsv.gz \
  --null_prob 0.001 \
  --min_coverage 5 \
  --fdr_threshold 0.05 \
  --threads 4
```

---

### 3.3 `prepare_matrix.sh`

Generates context-specific sliding-window matrices for two-group comparison. Takes per-sample `*_binomtest_result.tsv.gz` files from Sections 3.1–3.2 and integrates them across all samples and replicates within each window.

**Options:**

| Option | Default | Description |
|---|---|---|
| `--fasta` | required | Reference FASTA or FAI |
| `--group1` | required | Group 1 BinomTest TSV.gz files (multiple allowed) |
| `--group2` | required | Group 2 BinomTest TSV.gz files (multiple allowed) |
| `--group_labels` | `group1 group2` | Group names used in downstream analyses |
| `--window` | 300 | Window size (bp) |
| `--slide` | 200 | Slide size (bp) |
| `--output` | `./matrix_out` | Output directory |
| `--tmpdir` | system default | Temporary directory |

**Output format (`<g1>_<g2>_<ctx>_matrix.tsv.gz`):**

A 13-column layout derived from `bedtools intersect -wa -wb`, containing window coordinates, site coordinates, group/sample labels, strand, methylated counts, unmethylated counts, and coverage.

**Example (2 replicates per group):**

```bash
bash prepare_matrix.sh \
  --fasta /path/to/TAIR10.fasta \
  --group1 binom/WT_1_binomtest_result.tsv.gz \
           binom/WT_2_binomtest_result.tsv.gz \
  --group2 binom/MT_1_binomtest_result.tsv.gz \
           binom/MT_2_binomtest_result.tsv.gz \
  --group_labels WT MT \
  --window 500 \
  --slide 300 \
  --output prep_out
```

**Example (4 replicates per group):**

```bash
bash prepare_matrix.sh \
  --fasta /path/to/TAIR10.fasta \
  --group1 binom/WT_1_binomtest_result.tsv.gz \
           binom/WT_2_binomtest_result.tsv.gz \
           binom/WT_3_binomtest_result.tsv.gz \
           binom/WT_4_binomtest_result.tsv.gz \
  --group2 binom/MT_1_binomtest_result.tsv.gz \
           binom/MT_2_binomtest_result.tsv.gz \
           binom/MT_3_binomtest_result.tsv.gz \
           binom/MT_4_binomtest_result.tsv.gz \
  --group_labels WT MT \
  --window 500 \
  --slide 300 \
  --output prep_out
```

---

### 3.4 `run_glmmDMR.R`

Fits a GLMM to each window and estimates p-values, methylation differences (delta), and summary statistics. The framework supports four model configurations defined by the combination of response distribution (`--family`) and data representation (`--mode`).

**Model configurations:**

| `--family` | `--mode` | Response variable | Data representation |
|---|---|---|---|
| `beta` | `site` | Methylation proportion per cytosine | Site-level (**recommended**) |
| `beta` | `aggregate` | Window-mean methylation proportion | Aggregated |
| `binom` | `site` | Methylated read counts per cytosine | Site-level |
| `binom` | `aggregate` | Total methylated counts per window | Aggregated |

> **Recommendation:** `--family beta --mode site` is recommended based on benchmarking against simulated datasets with known ground-truth DMRs. This configuration maintains the highest precision across a broad range of methylation effect sizes, including subtle differences, and most effectively suppresses false positives driven by high replicate-level methylation variability. The beta distribution appropriately models the bounded [0,1] nature of methylation proportions, while site-level representation preserves cytosine-resolution information within each window.

In all configurations, group differences are modeled as fixed effects and biological replicates as random intercepts ($u_i \sim N(0, \sigma^2)$) when `--random_effect` is enabled. Statistical significance is assessed by likelihood ratio test comparing the full model against a reduced model without the group fixed effect, with p-values adjusted using the Benjamini–Hochberg method.

**Options:**

*Input / Output:*

| Option | Default | Description |
|---|---|---|
| `-i, --infile` | required | Input window matrix |
| `-o, --out_prefix` | required | Output prefix |

*Data selection:*

| Option | Default | Description |
|---|---|---|
| `-c, --context` | NULL | Context filter (`CpG` / `CHG` / `CHH`) |
| `--group1`, `--group2` | required | Comparison group labels |

*Model:*

| Option | Default | Description |
|---|---|---|
| `--family` | `beta` | Response distribution: `binom` or `beta` |
| `--mode` | `site` | Data representation: `aggregate` or `site` |
| `--random_effect` | TRUE | Model biological replicates as random intercepts |

*Pre-filter (applied before GLMM fitting):*

| Option | Default | Description |
|---|---|---|
| `--min_reps_g1`, `--min_reps_g2` | 2 | Minimum replicates per group required per window |
| `--min_cov` | 0 | Minimum coverage per cytosine site |
| `--min_sites_win` | 0 | Minimum cytosine sites per window |
| `--prefilter_delta` | 0 | Pre-remove windows with |delta| below this value (for speed) |

*Computation:*

| Option | Default | Description |
|---|---|---|
| `--workers` | 4 | Number of parallel workers |
| `--batches` | 50 | Number of batches for parallel processing |
| `--max_globals_mb` | 1000 | future globals size limit (MB) |
| `--seed` | 1 | Random seed for reproducibility |

**Output format (`*_fit_<family>_<mode>.tsv.gz`):**

| Column | Description |
|---|---|
| `chr` | Chromosome name |
| `start` | Window start position |
| `end` | Window end position |
| `model` | Fitted model identifier |
| `p` | Window-level p-value (BH-adjusted) |
| `delta` | Estimated methylation difference (group2 − group1) |
| `mean_rate1` | Mean methylation rate in group 1 |
| `mean_rate2` | Mean methylation rate in group 2 |
| `aic_diff` | AIC difference (full vs. reduced model) |
| `bic_diff` | BIC difference (full vs. reduced model) |

**Example:**

```bash
Rscript run_glmmDMR.R \
  -i prep_out/WT_MT_CpG_matrix.tsv.gz \
  -o glmm_out/WT_MT_CpG \
  --group1 WT --group2 MT \
  --family beta \
  --mode site \
  --min_reps_g1 2 --min_reps_g2 2 \
  --min_sites_win 3 \
  --min_cov 5 \
  --random_effect \
  --workers 8 --batches 200 --max_globals_mb 2000
```

---

### 3.5 `merge_window.R`

Integrates window-level GLMM statistics into DMRs using a seed-based region construction strategy. Window-level p-values within candidate regions are combined using the **Stouffer method**: each p-value is converted to a z-score ($z_i = \Phi^{-1}(1 - p_i)$) and the combined p-value is computed as $p_{\text{combined}} = 2\Phi(-|\sum z_i / \sqrt{k}|)$, where $k$ is the number of windows. This approach weights each window equally and yields a combined p-value that reflects both the number and consistency of significant windows within a region.

All windows are first assigned a methylation direction based on the sign of the estimated group effect (`delta`): positive delta → hypermethylated; negative delta → hypomethylated. Only windows sharing the same direction are integrated into a common region.

**Supported merge modes:**

| Mode | Description |
|---|---|
| `single_seed` | Starts from a single highly significant window (seed) and extends conservatively in both directions. Best for sharp, localized signals. |
| `multi_seed` | Defines seeds as clusters of neighboring significant windows evaluated jointly by combined p-value. Best for broad regions with multiple moderate signals. |
| `hybrid_seed` | Applies `multi_seed` genome-wide first, then applies `single_seed` to regions not covered by any multi-seed. **Default and recommended** for most use cases. |

**Options:**

*Input / Output:*

| Option | Default | Description |
|---|---|---|
| `--windows` | required | GLMM window results (`*_fit_<family>_<mode>.tsv.gz`) |
| `--out-prefix` | `results/dmr` | Output file prefix |

*Mode:*

| Option | Default | Description |
|---|---|---|
| `--merge-mode` | `hybrid_seed` | Region construction strategy |

*Seed / Extension (core detection parameters):*

| Option | Default | Description |
|---|---|---|
| `--p-seed` | 0.05 | Significance threshold for seed window definition |
| `--p-extend` | 0.05 | Significance threshold for region extension |
| `--max-gap-bp` | 200 | Maximum gap (bp) between adjacent windows for merging |
| `--min-windows` | 1 | Minimum number of windows required to report a DMR |
| `--min-delta` | 0 | Minimum effect size threshold during extension |
| `--max-p-degradation` | 1.2 | Maximum allowed increase in combined p-value during extension (1.0 = no worsening allowed) |
| `--max-final-p` | 1.0 | Maximum combined p-value of a retained DMR |
| `--min-strong-windows` | 0.5 | Minimum fraction of windows with p ≤ `--p-seed` in the final DMR |

*Adaptive effect size filter (recommended for `multi_seed` and `hybrid_seed`):*

| Option | Default | Description |
|---|---|---|
| `--adaptive-delta` | FALSE | Enable two-stage detection with data-driven effect size thresholding |
| `--adaptive-delta-method` | `median_ratio` | Threshold method: `median_ratio`, `q50`, `q25`, `q10`, `mad` |
| `--adaptive-delta-ratio` | 0.6 | Scaling ratio used in `median_ratio` method |

> **How adaptive delta works:** A first pass detects candidate DMRs without any effect size filter. The median absolute Δmethylation across all candidate DMRs is then computed, and a minimum |Δmethylation| threshold (set to `--adaptive-delta-ratio` × median) is applied in a second pass, removing extensions driven by weak or inconsistent methylation changes. `q25` is a practical starting point for initial parameter tuning.

*Multi-seed specific:*

| Option | Default | Description |
|---|---|---|
| `--seed-min-windows` | 1 | Minimum number of seed windows required for `multi_seed` and `hybrid_seed` |

*Post-filter (consistency checks after detection):*

| Option | Default | Description |
|---|---|---|
| `--post-filter` | FALSE | Enable post-detection quality filtering |
| `--min-median-p` | 0.01 | Maximum median p-value within a retained DMR |
| `--min-consistent-frac` | 0.5 | Minimum fraction of windows with p ≤ `--p-seed` |

*Length / median-p filters:*

| Option | Default | Description |
|---|---|---|
| `--min-dmr-length` | 0 | Minimum final DMR length (bp) |
| `--max-median-p` | 1.0 | Independent filter on DMR median p-value |

*Overlap merge:*

| Option | Default | Description |
|---|---|---|
| `--merge-overlaps` | FALSE | Re-merge overlapping/nearby same-direction DMRs after detection |
| `--merge-overlaps-gap` | 0 | Maximum gap (bp) allowed between DMRs during re-merge |

**Output format:**

| File | Pattern | Description |
|---|---|---|
| TSV | `*_dmrs_<mode>.tsv` | Full DMR table |
| BED | `*_dmrs_<mode>.bed` | Genome browser-compatible intervals |

**Main TSV output columns:**

| Column | Description |
|---|---|
| `chr` | Chromosome name |
| `start` | DMR start position |
| `end` | DMR end position |
| `n_windows` | Number of windows included in the DMR |
| `direction` | Direction of methylation change (`hyper` or `hypo`) |
| `combined_p` | Stouffer-combined p-value for the DMR |
| `delta` | Mean methylation difference across windows in the DMR |

**Example:**

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

### 3.6 `make_binned_methylation_bigwig.R` (optional)

Converts per-site methylation data from BinomTest output into bin-wise mean methylation bigWig files for genome browser visualization.

**Options:**

| Option | Default | Description |
|---|---|---|
| `-i, --input` | required | BinomTest result (`*_binomtest_result.tsv.gz`) |
| `-b, --binsize` | 50 | Bin size (bp) |
| `-o, --output` | required | Output bigWig file |
| `--genome` | required | Chromosome sizes file (chrom.sizes) |
| `--context` | `CpG` | Methylation context: `CpG`, `CHG`, or `CHH` |

**Example:**

```bash
Rscript make_binned_methylation_bigwig.R \
  -i binom/sample_binomtest_result.tsv.gz \
  -b 50 \
  -o wig/sample_CpG.bw \
  --genome /path/to/TAIR10.chromInfo \
  --context CpG
```

---

### 3.7 `make_binned_variance_bigwig.py` (optional)

Computes replicate-level methylation variance from bin-averaged values across multiple bigWig tracks and writes a variance bigWig file. This is useful for visualizing regions of high replicate-level variability alongside DMR calls.

**Options:**

| Option | Default | Description |
|---|---|---|
| `--inputs` | required | Input bigWig files (one per replicate) |
| `--output` | required | Output variance bigWig |
| `--bin-size` | 200 | Bin size (bp) |
| `--min-tracks` | 2 | Minimum number of tracks required per bin |
| `--norm` | `none` | Normalization: `none` or `log2p1` |

**Output:** Per-bin variance across input tracks (score in bigWig format).

**Example:**

```bash
python make_binned_variance_bigwig.py \
  --inputs wig/rep1_CpG.bw \
           wig/rep2_CpG.bw \
           wig/rep3_CpG.bw \
           wig/rep4_CpG.bw \
  --output wig/group_CpG.variance.bw \
  --bin-size 200 \
  --min-tracks 2 \
  --norm none
```

---

## 4. Troubleshooting

**1) No windows remain after filtering**

- Cause: Thresholds such as `--min_cov`, `--min_sites_win`, or `--min_reps_g1`/`--min_reps_g2` may be too strict for the dataset.
- Action: Relax thresholds stepwise. Check coverage distribution in your data before setting `--min_cov`.

**2) GLMM fitting is slow or runs out of memory**

- Action: Reduce `--workers`, increase `--batches`, and set `TMPDIR` to a fast storage location with sufficient capacity:

```bash
export TMPDIR=/path/to/large_storage/tmp
mkdir -p "$TMPDIR"
```

**3) Beta regression fails to converge for many windows**

- Cause: Beta regression can fail when methylation proportions are near 0 or 1, or when within-group variance is very low.
- Action: This is expected behavior; such windows are automatically excluded from output. If a large fraction of windows fail, consider switching to `--family binom` as a fallback, or check that `--min_sites_win` and `--min_cov` are set appropriately.

**4) Chromosome name mismatch during bigWig generation**

- Cause: Chromosome names in input files and the chrom.sizes file are inconsistent (e.g., `chr1` vs. `1`).
- Action: Harmonize chromosome naming across all input files and the reference chrom.sizes, then rerun.

**5) Too many or too few DMRs detected**

- Too many DMRs: Consider enabling `--adaptive-delta` or increasing `--max-final-p` stringency (lower value). Applying `--post-filter` with `--min-consistent-frac` can also reduce spurious detections.
- Too few DMRs: Relax `--max-final-p`, reduce `--p-seed`, or switch from `single_seed` to `multi_seed` or `hybrid_seed`.

**6) DMRs are highly fragmented**

- Cause: The `single_seed` strategy with strict parameters may produce many small disconnected regions.
- Action: Switch to `multi_seed` or `hybrid_seed`, increase `--max-gap-bp`, or enable `--merge-overlaps`.
