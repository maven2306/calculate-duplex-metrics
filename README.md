# Calculate Duplex Metrics

An R CLI tool to take in summarised read information and output a variety of duplex metrics.

### Available metrics

Individual metrics (selectable individually)
- frac_singletons
- efficiency
- drop_out_rate

Grouped metrics

GC metrics 
- gc_single
- gc_both
- gc_deviation

Family stats
- total_families
- family_mean
- family_median
- family_max
- families_gt1
- single_families
- paired_families
- paired_and_gt1

## Implementation overview

- The CLI entrypoint is `main.R`
- Argument parsing and validation are handled in `cli.R`
- Metric execution logic is in `calculate.R`
- Core metric implementations are defined in `R/calculate_nanoseq_functions.R`

Metric selection is resolved **before computation**.  
Only the requested individual metrics and/or metric groups are evaluated.

### GC metric behaviour
- GC metrics are computed only when a reference genome object (.fasta) is provided.
- `--metrics` defaults to `all`.
  - If `--ref_fasta` is not provided, GC metrics are skipped and a message is printed to the console.
  - If `--ref_fasta` is provided, GC metrics are computed (may return NA if insufficient data).
- If GC metrics are explicitly requested (e.g. `--metrics gc`) but no reference FASTA is supplied, the program exits with an error.


## Installation and Usage

### Option A: Using Docker (Recommended)

This method packages the script and all its dependencies into a self-contained environment. It is the most reliable way to run the analysis, as it guarantees that the exact same software versions are used every time.

#### Requirements
- **Docker:** You must have Docker installed and the Docker daemon running. You can download it from the [Docker website](https://www.docker.com/products/docker-desktop/).

#### Installation Steps

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/WEHIGenomicsRnD/calculate-duplex-metrics.git
    ```

2.  **Navigate to the project directory:**
    ```bash
    cd calculate-duplex-metrics
    ```

3.  **Build the Docker image:** This command reads the `Dockerfile` and builds a container image named `calculate-duplex-metrics`. This may take several minutes the first time you run it.
    ```bash
    docker build -t calculate-duplex-metrics .
    ```

4.  **(Optional) Verify the image:** You can check that the image was built successfully by listing your local Docker images.
    ```bash
    docker images
    ```
    You should see `calculate-duplex-metrics` in the list.

#### Default Usage Example

To run the tool, you use the `docker run` command. The `-v` flags are essential for allowing the Docker container to access files on your local machine.

```bash
docker run --rm \
  -v "$(pwd)/data:/app/data" \
  -v "$(pwd)/out:/app/out" \
  calculate-duplex-metrics \
  --input data/test.rinfo \
  --output out/default.csv
```

-   `-v "$(pwd)/data:/app/data"`: This "mounts" your local `data` directory into the `/app/data` directory inside the container, so the script can find the input file.
-   `-v "$(pwd)/out:/app/out"`: This mounts your local `out` directory into the `/app/out` directory inside the container, so the script can write the output file back to your machine.

### Option B: Local Installation (via `devtools`)

This method installs the required R packages directly onto your system. It is more flexible but less reproducible than Docker, as it will use your system's installed version of R and the latest available package versions.

#### Requirements
- **R:** R version 4.4.1 is recommended.
- **`devtools` R package:** This is used to install packages from GitHub.

#### Installation Steps

1.  **Install `devtools`:** If you don't have it already, open an R console and run:
    ```r
    install.packages("devtools")
    ```

2.  **Install Bioconductor dependencies:** This tool requires several packages from the Bioconductor repository.
    ```r
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install(c("Rsamtools", "GenomicRanges", "IRanges", "Biostrings"))
    ```

3.  **Install the package from GitHub:**
    ```r
    devtools::install_github('WEHIGenomicsRnD/calculate-duplex-metrics')
    ```

#### Default Usage Example

After installing the dependencies, you can run the script directly from your terminal within the cloned repository directory.

```bash
Rscript main.R \
  --input data/test.rinfo \
  --output out/default.csv
```

## Additional Usage Examples
#### Example: default mode with GC enabled (requires reference genome)

Note: The reference genome FASTA is user-provided and not included
in this repository. Any compatible reference genome may be used.

``` bash
Rscript main.R \
  --input data/NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001.txt \
  --output out/default_with_gc.csv \
  --ref_fasta ref/Escherichia_coli_ATCC_10798.fasta

```

#### Example: select individual metrics only

``` bash
Rscript main.R \
  --input data/NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001.txt \
  --output out/test_selected_metrics.csv \
  --metrics efficiency,drop_out_rate
```
Note: when listing multiple metrics, either omit spaces (efficiency,drop_out_rate) or quote the argument ("efficiency, drop_out_rate").

#### Example: select metric groups

``` bash
Rscript main.R \
  --input data/NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001.txt \
  --output out/test_family_metrics.csv \
  --metrics family
```

#### Example: mixed selection (individual + group)

``` bash
Rscript main.R \
  --input data/NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001.txt \
  --output out/test_mixed_metrics.csv \
  --metrics efficiency,family
```

#### Example: multiple input files

``` bash
Rscript main.R \
  --input data/a.txt data/b.txt \
  --output out/all_samples_metrics.csv \
  --metrics family \
  --cores 2

```


### CLI flags

``` bash
Required:
  -i, --input        One or more input rinfo files OR a directory containing rinfo files (.txt or .txt.gz)
                     Note: when --input is a directory, the tool selects matching files using --pattern (default: \.txt(\.gz)?$);
                           when --input is a list of files, --pattern is ignored
  -o, --output       Output CSV path (long format)

Optional:
  -s, --sample       Optional sample name(s). For multiple input files, provide
                     comma-separated names matching the number of files.
                     Note: if --input is a directory, --sample is not allowed.

      --pattern      Regex pattern used to select files when --input is a directory
                     (default: \.txt(\.gz)?$)

      --rlen         Read length (default: 151)
      --skips        Trimmed / ignored bases per read (NanoSeq = 5, xGen = 8)

      --ref_fasta    Reference genome FASTA (required for GC metrics)

      --metrics      Comma-separated list of metrics and/or metric groups
                     - Individual: frac_singletons, efficiency, drop_out_rate
                     - Groups: gc, family
                     (default: all)

      --cores        Number of CPU cores for parallel processing (default: 1)


```
Note: when listing multiple metrics, either omit spaces (efficiency,drop_out_rate) or quote the argument ("efficiency, drop_out_rate").


## Outputs

Output is written in long format:
```
sample,metric,value

``` 

#### Example:

```
sample,metric,value
NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001,frac_singletons,0.0418706803079419
NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001,efficiency,0.0490258329591602
NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001,drop_out_rate,0.320805646128878
NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001,total_families,23825702
NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001,family_mean,6.748161712309
NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001,family_median,5
NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001,family_max,50
NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001,families_gt1,16771629
NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001,single_families,6731955
NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001,paired_families,9994045
NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001,paired_and_gt1,8152302

```


#### Sanity check the CLI
```bash
Rscript main.R --help

```

## Testing 
Test that functions return valid numeric values, correct handling of edge cases (NA, zero reads, invalid inputs) and presence of expected metrics names.

#### Requirements
Packages: testthat

#### To run all tests
From the project root run:
```bash
Rscript tests/testthat.R

```
