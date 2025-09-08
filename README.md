# Calculate Duplex Metrics

An R CLI tool to take in summarised read information and output a variety of
QC duplex metrics.

**MVP metrics:** `efficiency`, `drop_out_rate`.

> Internals: the CLI delegates to `R/calc_duplex_metrics.R`, which sources
> `efficiency_nanoseq_functions.R` from the `G000204_duplex` repo.
> GC metrics are disabled for MVP to avoid extra genome/Bioc dependencies.

## Installation

## Requirements
- R ≥ 4.4 (lockfile was created on 4.5.1; both work)
- Packages: `argparse`, `magrittr`, `data.table`, `R.utils`
- Access to `G000204_duplex` repo (for `code/efficiency_nanoseq_functions.R`)

## Using renv (recommended)
```r
install.packages("renv", repos = "https://cloud.r-project.org")
renv::restore()
```
## Or install packages directly
```
install.packages(c("argparse","magrittr","data.table","R.utils"),
                 repos = "https://cloud.r-project.org")
```

## Usage

```bash
Rscript src/main.R
  --input  /abs/path/QC/read_info/<sample>.txt.gz \
  --output output/<sample>_duplex_metrics.csv \
  --rfunc_dir /abs/path/G000204_duplex/code/efficiency_nanoseq_functions.R \
  --rlen 151 --skips 5
```
### CLI flags
-i, --input rinfo file (.txt or .txt.gz) |
-o, --output output CSV path (long format) |
-s, --sample sample ID (defaults to input basename) |
--rfunc_dir folder OR file for efficiency_nanoseq_functions.R |
--rlen read length (default: 151) |
--skips trimmed/ignored bases per read (Nano=5, xGEN=8) |
-v, --verbose verbose logging

## Outputs
```
sample,metric,value
<sample>,efficiency,<0–1>
<sample>,drop_out_rate,<0–1>
```
