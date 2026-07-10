#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

show_help() {
    echo "Usage: prepare_matrix_group_based.sh --fasta <genome.fa or .fai> \\"
    echo "                             --group1 <file1.tsv.gz> <file2.tsv.gz> ... \\"
    echo "                             --group2 <file1.tsv.gz> <file2.tsv.gz> ... \\"
    echo "                             [--group_labels group1_label group2_label] \\"
    echo "                             [--window 1000] \\"
    echo "                             [--slide 500] \\"
    echo "                             [--output outdir] \\"
    echo "                             [--tmpdir tmpdir] \\"
    echo "                             [--threads N] \\"
    echo "                             [--help]"
    echo ""
    echo "Options:"
    echo "  --fasta            Fasta or .fai file used to generate genome windows"
    echo "  --group1           List of binomtest_result.tsv.gz files for group 1"
    echo "  --group2           List of binomtest_result.tsv.gz files for group 2"
    echo "  --group_labels     Labels for group1 and group2 (default: group1 group2)"
    echo "  --window           Sliding window size in bp (default: 1000)"
    echo "  --slide            Sliding window step in bp (default: 500)"
    echo "  --output           Output directory (default: ./matrix_out)"
    echo "  --tmpdir           Temporary directory for intermediate files"
    echo "  --threads          Threads for parallel decompress/sort/compress (default: 12)"
    echo "  --help             Show this message and exit"
    exit 0
}

# Default values
window=300
slide=200
outdir="./matrix_out"
tmpdir=""
group1_label="group1"
group2_label="group2"
threads=12
group1_files=()
group2_files=()

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --fasta) fasta="$2"; shift 2 ;;
        --group1) shift; while [[ $# -gt 0 && "$1" != --* ]]; do group1_files+=("$1"); shift; done ;;
        --group2) shift; while [[ $# -gt 0 && "$1" != --* ]]; do group2_files+=("$1"); shift; done ;;
        --group_labels) group1_label="$2"; group2_label="$3"; shift 3 ;;
        --window) window="$2"; shift 2 ;;
        --slide) slide="$2"; shift 2 ;;
        --output) outdir="$2"; shift 2 ;;
        --tmpdir) tmpdir="$2"; shift 2 ;;
        --threads) threads="$2"; shift 2 ;;
        --help) show_help ;;
        *) echo "Unknown option $1"; exit 1 ;;
    esac
done

# Check required options
if [[ -z "$fasta" || ${#group1_files[@]} -eq 0 || ${#group2_files[@]} -eq 0 ]]; then
    echo "Error: --fasta, --group1, and --group2 are required."
    echo "Use --help to see usage."
    exit 1
fi

# Validate input files
for f in "${group1_files[@]}" "${group2_files[@]}"; do
    [[ -f "$f" ]] || { echo "Error: input file not found: $f"; exit 1; }
done

# Compute common prefix/suffix across all sample basenames
all_bases=()
for f in "${group1_files[@]}" "${group2_files[@]}"; do
    all_bases+=("$(basename "$f" .tsv.gz)")
done

common_prefix() {
    local arr=("$@")
    local prefix="${arr[0]}"
    for s in "${arr[@]}"; do
        while [[ -n "$prefix" && "${s#"$prefix"}" == "$s" ]]; do
            prefix="${prefix%?}"
        done
    done
    echo "$prefix"
}

common_suffix() {
    local arr=("$@")
    local suffix="${arr[0]}"
    for s in "${arr[@]}"; do
        while [[ -n "$suffix" && "${s%"$suffix"}" == "$s" ]]; do
            suffix="${suffix#?}"
        done
    done
    echo "$suffix"
}

common_pref=$(common_prefix "${all_bases[@]}")
common_suf=$(common_suffix "${all_bases[@]}")
echo "[INFO] Common sample prefix: '${common_pref}'"
echo "[INFO] Common sample suffix: '${common_suf}'"

# Prepare output and temp dirs
mkdir -p "$outdir"
if [[ -z "$tmpdir" ]]; then
    tmpdir=$(mktemp -d)
    trap 'rm -rf "$tmpdir"' EXIT
else
    mkdir -p "$tmpdir"
fi

# Check required tools
command -v bedtools >/dev/null 2>&1 || { echo "Error: bedtools not found"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "Error: samtools not found"; exit 1; }

# Prefer pigz for parallel (de)compression; fall back to gzip/zcat if unavailable.
if command -v unpigz >/dev/null 2>&1; then
    DECOMP=(unpigz -c)
elif command -v zcat >/dev/null 2>&1; then
    DECOMP=(zcat)
else
    echo "Error: neither unpigz nor zcat found"; exit 1
fi
have_pigz=0
command -v pigz >/dev/null 2>&1 && have_pigz=1

# Per-context share of threads (the 3 contexts are processed concurrently).
ctx_threads=$(( threads / 3 )); (( ctx_threads < 1 )) && ctx_threads=1
sort_buf="2G"
echo "[INFO] threads=$threads (per-context sort/compress threads=$ctx_threads), decomp='${DECOMP[*]}', pigz=$have_pigz"

# Prepare .fai if needed
if [[ "$fasta" == *.fai ]]; then
    fai="$fasta"
else
    fai="${fasta}.fai"
    [[ -f "$fai" ]] || samtools faidx "$fasta"
fi

# Create window bed file
window_bed="$tmpdir/windows.bed"
bedtools makewindows -g "$fai" -w "$window" -s "$slide" > "$window_bed"
if [[ ! -s "$window_bed" ]]; then
    echo "Error: window bed is empty. Check fasta/fai and window/slide settings."
    exit 1
fi

# --- Phase 1: one decompress pass per sample, split into per-context BEDs ---
# Each input file is read exactly once and demultiplexed by context ($6) into
# the CpG / CHG / CHH BEDs, instead of re-reading the whole file once per
# context. Samples are independent, so they run in parallel.
process_sample() {
    local label="$1" f="$2"
    local base sample_id
    base=$(basename "$f" .tsv.gz)
    sample_id="$base"
    [[ -n "$common_pref" ]] && sample_id="${sample_id#"$common_pref"}"
    [[ -n "$common_suf" ]] && sample_id="${sample_id%"$common_suf"}"
    [[ -z "$sample_id" ]] && sample_id="$base"
    echo "[INFO] Splitting by context: $label / $sample_id"
    "${DECOMP[@]}" "$f" | awk -v OFS='\t' -v g="$label" -v s="$sample_id" \
        -v f_cpg="$tmpdir/${label}_${base}_CpG_cytosines.bed" \
        -v f_chg="$tmpdir/${label}_${base}_CHG_cytosines.bed" \
        -v f_chh="$tmpdir/${label}_${base}_CHH_cytosines.bed" '
        NR > 1 {
            cov = $4 + $5
            # pos ($2) is 1-based (Bismark); BED is 0-based half-open -> [pos-1, pos)
            line = $1 OFS ($2 - 1) OFS $2 OFS $6 OFS g OFS s OFS $3 OFS $4 OFS $5 OFS cov
            if      ($6 == "CpG") print line > f_cpg
            else if ($6 == "CHG") print line > f_chg
            else if ($6 == "CHH") print line > f_chh
        }
    '
}

sample_entries=()
for f in "${group1_files[@]}"; do sample_entries+=("$group1_label"$'\t'"$f"); done
for f in "${group2_files[@]}"; do sample_entries+=("$group2_label"$'\t'"$f"); done

echo "[INFO] Phase 1: splitting ${#sample_entries[@]} samples by context (up to $threads in parallel)"
n_samples=${#sample_entries[@]}
fail=0
for (( i = 0; i < n_samples; i += threads )); do
    pids=()
    for (( j = i; j < i + threads && j < n_samples; j++ )); do
        process_sample "${sample_entries[j]%%$'\t'*}" "${sample_entries[j]#*$'\t'}" &
        pids+=($!)
    done
    for p in "${pids[@]}"; do wait "$p" || fail=1; done
done
(( fail == 0 )) || { echo "Error: a sample-splitting job failed"; exit 1; }

# --- Phase 2: per-context sort + window intersect (contexts run in parallel) ---
# The sites are sorted into .fai chromosome order (external merge sort, spills to
# disk, so memory stays bounded even at ~1e8 rows) so that
# 'bedtools intersect -sorted -g' streams both inputs (chromsweep) instead of
# loading -b into an in-memory interval tree. A leading column with each
# chromosome's fai rank makes the numeric sort reproduce fai order; it is
# stripped again before intersect.
process_context_intersect() {
    local ctx="$1"
    local -a beds comp
    mapfile -t beds < <(find "$tmpdir" -maxdepth 1 -name "*_${ctx}_cytosines.bed" -print)
    if (( ${#beds[@]} == 0 )); then
        echo "[WARNING] No $ctx cytosine BED files found"
        return 0
    fi
    if (( have_pigz )); then comp=(pigz -p "$ctx_threads"); else comp=(gzip); fi
    echo "[INFO] Phase 2: $ctx — sorting ${#beds[@]} BEDs (fai order) + window intersect"
    cat "${beds[@]}" \
      | awk 'NR==FNR { rank[$1] = FNR; next } { print rank[$1] "\t" $0 }' "$fai" - \
      | LC_ALL=C sort -k1,1n -k3,3n -S "$sort_buf" --parallel="$ctx_threads" -T "$tmpdir" \
      | cut -f2- \
      | bedtools intersect -sorted -g "$fai" -a "$window_bed" -b - -wa -wb \
      | "${comp[@]}" > "$outdir/${group1_label}_${group2_label}_${ctx}_matrix.tsv.gz"
    rm -f "${beds[@]}"
}

echo "[INFO] Phase 2: intersecting 3 contexts in parallel"
pids=()
for ctx in CpG CHG CHH; do
    process_context_intersect "$ctx" &
    pids+=($!)
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done
(( fail == 0 )) || { echo "Error: a context intersect job failed"; exit 1; }

echo "[DONE] Matrix files written to $outdir"
