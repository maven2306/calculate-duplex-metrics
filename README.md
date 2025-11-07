# Calculate Duplex Metrics

An R CLI tool to take in summarised read information and output a variety of duplex metrics.

**QC metrics:** `efficiency`, `drop_out_rate`.

> Internals: the CLI entrypoint is `src/main.R` (argument parsing + I/O). It delegates computation to `R/calc_duplex_metrics.R`, which sources from the `src/efficiency_nanoseq_functions.R`. GC metrics are disabled for MVP to avoid extra genome/Bioc dependencies.

## Installation

-   R version 4.4.1

-   Packages: `argparse`, `magrittr`, `data.table`, `R.utils`

## Requirement

-   Access to `G000204_duplex` repo (for `code/efficiency_nanoseq_functions.R`) then move `efficiency_nanoseq_functions.R` to `src` within calculate-duplex-metrics repo.

## Using renv (recommended)

```         
R -q -e "install.packages('renv', repos = '<https://cloud.r-project.org>'); renv::restore()"
```

## Or install directly (not use renv)

```         
install.packages(c("argparse","magrittr","data.table","R.utils"),
                 repos = "https://cloud.r-project.org")
```

## Usage

#### Example:

``` bash
Rscript src/main.R \
  --input  "/Users/dylandinh/Desktop/QC/read_info/NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001.txt.gz" \
  --output "/Users/dylandinh/calculate-duplex-metrics/output/NanoMB1Rep1_duplex_metrics.csv" \
  --sample "NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001" \
  --rfunc_dir "/Users/dylandinh/calculate-duplex-metrics/src/efficiency_nanoseq_functions.R" \
  --rlen 151 --skips 5
```

#### Minimal version:

``` bash
Rscript src/main.R
  --input  /abs/path/QC/read_info/<sample>.txt.gz \
  --output output/<sample>_duplex_metrics.csv \
  --rfunc_dir /abs/path/calculate-duplex-metrics/src/efficiency_nanoseq_functions.R \
  --rlen 151 --skips 5
```

### CLI flags

-i, --input : rinfo file (.txt or .txt.gz) | -o, --output : output CSV path (long format) | -s, --sample : sample ID (defaults to input basename) | --rfunc_dir : folder OR file for efficiency_nanoseq_functions.R | --rlen : read length (default: 151) | --skips : trimmed/ignored bases per read (Nano=5, xGEN=8) | -v, --verbose : verbose logging

#### Sanity check the CLI
```bash
Rscript src/main.R --help
```
## Outputs

```         
sample,metric,value
<sample>,efficiency,<0–1>
<sample>,drop_out_rate,<0–1>
```

#### Example:

```         
sample,metric,value
NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001,efficiency,0.0490258329591602
NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001,drop_out_rate,0.320805646128878
```
