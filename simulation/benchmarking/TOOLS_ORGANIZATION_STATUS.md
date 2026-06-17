# ベンチマーキングツール整理状況

## 現在の状態 (2026-06-17)

### ✅ GitHub に個別スクリプトがある (2/6)

| ツール | ファイル | 状態 |
|--------|---------|------|
| DSS | `04.run_DSS.R` | ✅ 統合済 |
| methylKit | `04.run_methylKit.R` | ✅ 統合済 |

### ⚠️ リモートサーバーにのみ存在 (4/6)

| ツール | 場所 | 実装方法 |
|--------|------|---------|
| Fisher | `bash_variance.sh` (行285-306) | bash 統合版 |
| metilene | `bash_variance.sh` (行309-317) | bash 統合版 |
| DMRfinder | `bash_variance.sh` (行320-346) | bash 統合版 |
| MACAU2 | `bash_variance.sh` (行348-359) | bash 統合版 |

---

## 問題点

**Fisher, metilene, DMRfinder, MACAU2 は個別の R スクリプトとして GitHub に存在しません**

- 実行ロジックはすべて `bash_variance.sh` に埋め込まれている
- 個別テスト・カスタマイズが困難
- 他のツール（DSS, methylKit）との統一性がない

---

## 推奨解決策

### オプション A: 個別 R スクリプトを作成 (推奨)

`bash_variance.sh` から以下を抽出して個別ファイルを作成：

#### 1. `04.run_fisher.R`
```bash
# 行287-291 の run_fisher.R コマンドを R スクリプトに
Input: window データ
Output: output_for_fisher/fisher_out_dmrs.tsv
```

#### 2. `04.run_metilene.R`
```bash
# 行311-314 の metilene コマンドを R ラッパーに
Input: sites_CG_forMetilene.txt
Output: output_for_metilene/metilene_out.tsv
```

#### 3. `04.run_dmrfinder.R`
```bash
# 行324-346 の DMRfinder Python+R コマンドをラッパーに
Input: sites_CG_forDMRfinder_*.txt
Output: output_for_DMRfinder/out_findDMRs.txt
```

#### 4. `04.run_macau2.R`
```bash
# 行355 の 04.run_MACAU2.R を output_for_MACAU/ から取得
Input: サイトレベルデータ
Output: output_for_MACAU/MACAU_dmrs.tsv
```

**利点:**
- ✅ DSS/methylKit との統一性
- ✅ 個別テスト可能
- ✅ パラメータ調整が容易
- ✅ GitHub 単独で実行可能

---

### オプション B: bash_variance.sh を GitHub に追加

`simulation/benchmarking/` に追加：
- `bash_variance.sh` — 全ツール統合スクリプト
- `bash_variance_moderate.sh`
- `bash_variance_extreme.sh`

**利点:**
- 既存コードをそのまま使用
- 実行は簡単

**欠点:**
- ✗ bash に依存
- ✗ ツール毎のカスタマイズ困難
- ✗ 4 つのツールがドキュメント化されていない

---

## リモートサーバーの bash_variance.sh 内容

**ツール実行箇所:**

```bash
# 行285-306: Fisher's exact test
/usr/bin/time -v -o time_fisher.txt run_fisher.R \
  -i windows_CG.tsv.gz \
  -o output_for_fisher/fisher_out.tsv \
  -a WT -b MT

# 行309-317: metilene
/usr/bin/time -v -o time_metilene.txt /home/epigenome/.local/atools/metilene_v0.2-9/metilene \
  -i ./output_for_metilene/sites_CG_forMetilene.txt > ./output_for_metilene/metilene_out.tsv

# 行320-346: DMRfinder  
python /home/epigenome/.local/atools/DMRfinder/combine_CpG_sites.py \
  -o results.csv \
  sites_CG_forDMRfinder_*.txt
Rscript /home/epigenome/.local/atools/DMRfinder/findDMRs_fixed.r \
  --input results.mod.csv --output out_findDMRs.txt

# 行348-359: MACAU 2.0
Rscript 04.run_MACAU2.R > log_macau_output.txt 2>&1
```

---

## 次のアクション

**推奨 (オプション A の実施):**

1. `bash_variance.sh` から各ツール実行コマンドを抽出
2. 個別 R スクリプトを作成 (またはリモートから取得)
3. 相対パスに統一
4. `benchmarking/README.md` を更新してすべてのツールをドキュメント化
5. テスト実行して動作確認

**現在の状態:**
- ✅ DSS, methylKit: 個別スクリプト化済
- ⏳ Fisher, metilene, DMRfinder, MACAU2: 要抽出

---

## ファイル構成（完成時）

```
simulation/benchmarking/
├── README.md (更新必要)
├── 03.convert_sites_for_otherSoft.R
├── 04.run_DSS.R
├── 04.run_methylKit.R
├── 04.run_fisher.R (NEW)
├── 04.run_metilene.R (NEW)
├── 04.run_dmrfinder.R (NEW)
├── 04.run_macau2.R (NEW)
└── 05.evaluate_dmrs.R
```

---

## 参考: リモートサーバー構成

- **bash_variance.sh**: マスター実行スクリプト (643 行)
  - シミュレーション、データ変換、全ツール実行を統合
  - 他バリアント: bash_variance_moderate.sh, bash_variance_extreme.sh

- **output_* ディレクトリ**: ツール別出力
  - output_for_DSS/
  - output_for_methylKit/
  - output_for_fisher/
  - output_for_metilene/
  - output_for_DMRfinder/
  - output_for_MACAU/
