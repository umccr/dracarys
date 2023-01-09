
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

### 🌈 multiqc

- Generate a ‘tidier’ form of the `multiqc_data.json`
  [MultiQC](https://multiqc.info/) data summary file. Builds on
  functionality from
  [TidyMultiqc](https://github.com/multimeric/TidyMultiqc). TSV and/or
  Parquet outputs are generated. For the TSV file, each row corresponds
  to a single sample, and each column corresponds to a single quality
  metric/variable. See the [CLI](#cli) section below for options.

### 🌶 tso

- Generate a ‘tidier’ form of the UMCCR TSO500 ctDNA workflow results.
  Input is a GDS (or local) directory with ‘raw’ TSO500 ctDNA results,
  and output includes the relevant GDS files synced, their tidy form in
  parquet/tsv format, a HTML report with those results in table form,
  and an RDS R file that feeds into the HTML report. See the [CLI](#cli)
  section below for options.

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
    usage: dracarys.R tso [-h] -i IN_DIR -o OUT_DIR [-r REPORT_DIR] -p PREFIX
                          [-t TOKEN] [-g GDS_LOCAL_DIR] [--rds_dir RDS_DIR]
                          [-f {tsv,parquet,both}] [-n] [-q] [--quiet_rmd]

    optional arguments:
      -h, --help            show this help message and exit
      -i IN_DIR, --in_dir IN_DIR
                            💩 Directory with TSO500 ctDNA workflow results (GDS or
                            local).
      -o OUT_DIR, --out_dir OUT_DIR
                            🎁 Output tidy results to this directory.
      -r REPORT_DIR, --report_dir REPORT_DIR
                            ✨ Output HTML report with tidy RDS object to this
                            directory.
      -p PREFIX, --prefix PREFIX
                            💃 Prefix string (used for all results).
      -t TOKEN, --token TOKEN
                            🙈 ICA access token (def. ICA_ACCESS_TOKEN env var).
      -g GDS_LOCAL_DIR, --gds_local_dir GDS_LOCAL_DIR
                            📋 If input is a GDS directory, download the
                            'recognisable' files to this directory. If not
                            specified, files will be downloaded to
                            'out_dir/dracarys_gds_sync'.
      --rds_dir RDS_DIR     💧 Directory to save RDS object with results from tidy
                            function.
      -f {tsv,parquet,both}, --format {tsv,parquet,both}
                            🍦 Format of output (default: tsv).
      -n, --dryrun          🐫 Dry run (just print tibble with files to be tidied).
      -q, --quiet           😴 Shush all the logs (also see --quiet_rmd).
      --quiet_rmd           😴 Shush just the Rmd rendering logs.
