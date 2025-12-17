# Calculate Duplex Metrics

An R CLI tool to take in summarised read information and output a variety of duplex metrics.

### Available metrics

- frac_singletons
- efficiency
- drop_out_rate
- gc_single
- gc_both
- gc_deviation
- total_families
- family_mean
- family_median
- family_max
- families_gt1
- single_families
- paired_families
- paired_and_gt1


> Internals: the CLI entrypoint is `main.R`, which invokes `cli.R` for
> argument parsing and validation. Computation is orchestrated in
> `calculate.R` and delegated to `R/calc_duplex_metrics.R`, which sources
> metric logic from `R/efficiency_nanoseq_functions.R`.
>
> GC metrics are computed when a reference FASTA is provided; otherwise
> they are either set to `NA` (when GC is skipped) or the program errors
> if GC was requested but no reference was given.

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
  --input data/test.rinfo \
  --output test_duplex_metrics.csv \
  --sample test \
  --skip_gc TRUE

```

#### Example (with GC enabled, requires reference genome):

``` bash
Rscript main.R \
  --input data/test.rinfo \
  --output test_duplex_metrics_gc.csv \
  --sample test \
  --ref_fasta ref/Escherichia_coli_ATCC_10798.fasta \
  --skip_gc FALSE

```

#### Example (select specific metrics)

``` bash
Rscript main.R \
  --input data/test.rinfo \
  --output test_selected_metrics.csv \
  --sample test \
  --metrics efficiency,drop_out_rate \
  --skip_gc TRUE
```

#### Example (multiple input files)

``` bash
Rscript main.R \
  --input data/a.rinfo data/b.rinfo \
  --output all_samples_metrics.csv \
  --skip_gc TRUE \
  --cores 2
```

#### Example (input directory + pattern)
``` bash
Rscript main.R \
  --input_dir data \
  --pattern "\\.rinfo(\\.gz)?$" \
  --output all_samples_metrics.csv \
  --skip_gc TRUE

```

### CLI flags

``` bash
Required:
  -i, --input        One or more input rinfo files (.txt or .txt.gz)
      --input_dir    Directory containing rinfo files
  -o, --output       Output CSV path (long format, MultiQC-compatible) 
  
Optional:
  -s, --sample       Sample ID (only valid for a single input file, otherwise sample names are derived from filenames)
      --pattern      Regex pattern used with --input_dir
      --rlen         Read length (default: 151)
      --skips        Trimmed / ignored bases per read (NanoSeq = 5, xGen = 8)

      --ref_fasta    Reference genome FASTA (enables GC metrics)
      --skip_gc      TRUE / FALSE; disable GC even if FASTA is provided
                     (default: TRUE)

      --metrics      Optional comma-separated list of metrics to output
                     (default: all metrics returned by calculate_metrics_single)
      --cores        Number of CPU cores for parallel processing
                     (default: 1)
                     
  -v, --verbose      Verbose output

  
```
Note: Specify either --input or --input_dir (not both). When multiple input files are provided, results are combined into one output CSV.


#### Sanity check the CLI
```bash
Rscript main.R --help

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
