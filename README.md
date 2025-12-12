# Calculate Duplex Metrics

An R CLI tool to take in summarised read information and output a variety of duplex metrics.

**QC metrics now supported:** `frac_singletons`, `efficiency`, `drop_out_rate`, GC metrics, and family statistics.

> Internals: the CLI entrypoint is `main.R` (argument parsing + I/O). It delegates computation to `R/calc_duplex_metrics.R`, which sources from `R/efficiency_nanoseq_functions.R`.  
> GC metrics are computed when a reference FASTA is provided; otherwise they are either set to `NA` (when GC is skipped) or the script errors if GC was requested but no reference was given.

## Installation

- R version 4.4.1
- Packages: `argparse`, `magrittr`, `data.table`, `R.utils`, `Biostrings`,  
  `GenomicRanges`, `IRanges`, `Rsamtools`

## Using renv (recommended)

```r
R -q -e "install.packages('renv', repos = 'https://cloud.r-project.org'); renv::restore()"

```  

## Or install manually (not use renv)

```r
install.packages(
  c('argparse', 'magrittr', 'data.table', 'R.utils'),
  repos = 'https://cloud.r-project.org'
)

if (!requireNamespace('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager', repos = 'https://cloud.r-project.org')
}
BiocManager::install(c('Biostrings', 'GenomicRanges', 'IRanges', 'Rsamtools'))
```  

## Usage

#### Example (no GC computation):

``` bash
Rscript main.R \
  --input  data/test.rinfo \
  --output test_duplex_metrics.csv \
  --sample test \
  --rfunc_dir R/efficiency_nanoseq_functions.R \
  --rlen 151 \
  --skips 5 \
  --skip_gc TRUE
```

#### Example (with GC enabled, requires refrence genome):

``` bash
Rscript main.R \
  --input data/test.rinfo \
  --output test_duplex_metrics.csv \
  --sample test \
  --rfunc_dir R/efficiency_nanoseq_functions.R \
  --rlen 151 \
  --skips 5 \
  --ref_fasta ref/Ecoli.fa \
  --skip_gc FALSE
```


### CLI flags

```
Options:
  -i, --input        rinfo file (.txt or .txt.gz)
  -o, --output       output CSV path (long format)
  -s, --sample       sample ID (defaults to input filename)
  --rfunc_dir        directory or file for efficiency_nanoseq_functions.R
  --rlen             read length (default: 151)
  --skips            number of trimmed/ignored bases (Nano=5, xGEN=8)
  --ref_fasta        optional reference genome FASTA (enables GC metrics)
  --skip_gc          TRUE/FALSE; force-disable GC even if FASTA is provided
  -v, --verbose      verbose logging output

```

#### Sanity check the CLI
```bash
Rscript src/main.R --help
```
## Outputs

```         
sample,metric,value
<sample_id>,frac_singletons,<num>
<sample_id>,efficiency,<num>
<sample_id>,drop_out_rate,<num>
<sample_id>,gc_single,<num_or_NA>
<sample_id>,gc_both,<num_or_NA>
<sample_id>,gc_deviation,<num_or_NA>
<sample_id>,total_families,<int>
<sample_id>,family_mean,<num>
<sample_id>,family_median,<num>
<sample_id>,family_max,<int>
<sample_id>,families_gt1,<int>
<sample_id>,single_families,<int>
<sample_id>,paired_families,<int>
<sample_id>,paired_and_gt1,<int>

```

#### Example:

```         
sample,metric,value
test,frac_singletons,0.089
test,efficiency,0.064
test,drop_out_rate,-0.033
test,gc_single,NA
test,gc_both,NA
test,gc_deviation,NA
test,total_families,9
test,family_mean,5
test,family_median,3
test,family_max,13
test,families_gt1,5
test,single_families,4
test,paired_families,4
test,paired_and_gt1,3

```
