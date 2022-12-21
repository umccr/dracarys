
<!-- README.md is generated from README.Rmd. Please edit that file -->

# 🔥 dracarys - DRAGEN Workflow Post-Processing

![](https://emojis.slackmojis.com/emojis/images/1643517245/32837/dragon_wiggle.gif?1643517245 "Dragon Dance")

- Docs: <https://umccr.github.io/dracarys/>

[![Conda
install](https://anaconda.org/umccr/r-dracarys/badges/latest_release_date.svg)](https://anaconda.org/umccr/r-dracarys)

## Installation

``` r
remotes::install_github("umccr/dracarys")
```

- Or if used inside a conda environment:

``` bash
conda install r-dracarys -c umccr -c conda-forge -c bioconda
```

## Main modules

### ✨ multiqc

- Generate a ‘tidier’ form of the `multiqc_data.json` MultiQC summary
  file. TSV and Parquet outputs are generated. For the TSV file, each
  row corresponds to a single sample, and each column corresponds to a
  single quality metric/variable. See the [CLI](#cli) section below for
  options.

### tso

- Generate a ‘tidier’ form of the TSO500 ctDNA workflow results.

## 💻 CLI

A `dracarys` command line interface is available for convenience.

- If you’re using the conda package, the `dracarys.R` command will
  already be set up inside an activated conda environment.
- If you’re *not* using the conda package, you need to export the
  `dracarys/inst/cli/` directory to your `PATH` in order to use
  `dracarys.R`.

``` bash
dracarys_cli=$(Rscript -e 'x = system.file("cli", package = "dracarys"); cat(x, "\n")' | xargs)
export PATH="${dracarys_cli}:${PATH}"
```

    $ dracarys.R --version
    dracarys.R 0.5.0

    #------- multiqc -------#


    $ dracarys.R multiqc --help
    usage: dracarys.R multiqc [-h] -j JSON -o OUTDIR -p PREFIX [-q]
                              [-f {tsv,parquet,both}]

    optional arguments:
      -h, --help            show this help message and exit
      -j JSON, --json JSON  💩 Path to 'multiqc_data.json'.
      -o OUTDIR, --outdir OUTDIR
                            🎁 Output results to this directory.
      -p PREFIX, --prefix PREFIX
                            ✨ Prefix output files with this string.
      -q, --quiet           😴 Shush all the logs.
      -f {tsv,parquet,both}, --format {tsv,parquet,both}
                            🍦 Format of output (default: tsv).
    #------- tso -------#


    $ dracarys.R tso --help
    usage: dracarys.R tso [-h] -i INDIR -o OUTDIR [-q] [-n] [-f {tsv,parquet,rds}]

    optional arguments:
      -h, --help            show this help message and exit
      -i INDIR, --indir INDIR
                            💩 Path to directory with TSO ctDNA workflow results.
      -o OUTDIR, --outdir OUTDIR
                            🎁 Output tidy results to this directory.
      -q, --quiet           😴 Shush all the logs.
      -n, --dryrun          🐫 Dry run.
      -f {tsv,parquet,rds}, --format {tsv,parquet,rds}
                            🍦 Format of output (default: tsv).
