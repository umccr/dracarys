
<!-- README.md is generated from README.Rmd. Please edit that file -->

# 🔥 dracarys - DRAGEN Workflow Post-Processing

![](https://emojis.slackmojis.com/emojis/images/1643517245/32837/dragon_wiggle.gif?1643517245 "Dragon Dance")

- Docs: <https://umccr.github.io/dracarys/>

[![Conda
install](https://anaconda.org/umccr/r-dracarys/badges/latest_release_date.svg)](https://anaconda.org/umccr/r-dracarys)

## Goal

Given a GDS (or local) directory with results from a UMCCR workflow,
{dracarys} will grab files of interest and transform them into ‘tidier’
structures for output into TSV/Parquet format. For supported files see
the [Supported Files](#supported-files) section below.

## Installation

``` r
remotes::install_github("umccr/dracarys")
```

- Or if used inside a conda environment:

``` bash
conda install r-dracarys -c umccr -c conda-forge -c bioconda
```

## Supported Files

### 🐉 DRAGEN

- Generates a ‘tidier’ form of the UMCCR DRAGEN Somatic workflow
  outputs:
  - `wgs_contig_mean_cov_<phenotype>.csv`
  - `wgs_coverage_metrics_<phenotype>.csv`
  - `wgs_fine_hist_<phenotype>.csv`
  - `fragment_length_hist.csv`
  - `mapping_metrics.csv`
  - `ploidy_estimation_metrics.csv`
  - `replay.json`
  - `time_metrics.csv`
  - `vc_metrics.csv`

### 🌶 tso

- Generates a ‘tidier’ form of the UMCCR TSO500 ctDNA workflow outputs:
  - `SampleAnalysisResults.json.gz`
  - `AlignCollapseFusionCaller_metrics.json.gz`
  - `TMB_Trace.tsv`
  - `Fusions.csv`
  - `tmb.json.gz`
  - `msi.json.gz`
  - `fragment_length_hist.json.gz`
  - `TargetRegionCoverage.json.gz`

### 🌈 MultiQC

- Generates a ‘tidier’ form of the `multiqc_data.json`
  [MultiQC](https://multiqc.info/) data summary file. Builds on
  functionality from
  [TidyMultiqc](https://github.com/multimeric/TidyMultiqc). TSV and/or
  Parquet outputs are generated. For the TSV file, each row corresponds
  to a single sample, and each column corresponds to a single quality
  metric/variable.

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
    dracarys.R 0.6.0

    #------- tidy -------#


    $ dracarys.R tidy --help
    usage: dracarys.R tidy [-h] -i IN_DIR -o OUT_DIR -p PREFIX [-t TOKEN]
                           [-g GDS_LOCAL_DIR] [-f {tsv,parquet,both}] [-n] [-q]

    optional arguments:
      -h, --help            show this help message and exit
      -i IN_DIR, --in_dir IN_DIR
                            ⛄️ Directory with untidy UMCCR workflow results (GDS
                            or local).
      -o OUT_DIR, --out_dir OUT_DIR
                            🔥 Directory to output tidy results.
      -p PREFIX, --prefix PREFIX
                            🎻 Prefix string (used for all results).
      -t TOKEN, --token TOKEN
                            🙈 ICA access token (def. ICA_ACCESS_TOKEN env var).
      -g GDS_LOCAL_DIR, --gds_local_dir GDS_LOCAL_DIR
                            📋 If input is a GDS directory, download the
                            'recognisable' files to this directory. If not
                            specified, files will be downloaded to
                            '<out_dir>/dracarys_gds_sync'.
      -f {tsv,parquet,both}, --format {tsv,parquet,both}
                            🍦 Format of output (default: tsv).
      -n, --dryrun          🐫 Dry run (just show files to be tidied).
      -q, --quiet           😴 Shush all the logs.
