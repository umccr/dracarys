
- <a href="#-dracarys---dragen-workflow-post-processing-"
  id="toc--dracarys---dragen-workflow-post-processing-">🐲 dracarys -
  DRAGEN Workflow Post-Processing 🔥</a>
  - <a href="#installation" id="toc-installation">Installation</a>
  - <a href="#main-modules" id="toc-main-modules">Main modules</a>
    - <a href="#id_-tidy" id="toc-id_-tidy">✨ Tidy</a>
  - <a href="#id_-cli" id="toc-id_-cli">💻 CLI</a>

<!-- README.md is generated from README.Rmd. Please edit that file -->

# 🐲 dracarys - DRAGEN Workflow Post-Processing 🔥

- Docs: <https://umccr.github.io/dracarys/>

[![Conda
install](https://anaconda.org/umccr/r-dracarys/badges/installer/conda.svg)](https://anaconda.org/umccr/r-dracarys)

## Installation

``` r
remotes::install_github("umccr/dracarys")
```

- Or if used inside a conda environment:

``` bash
conda install r-dracarys -c umccr -c conda-forge -c bioconda
```

## Main modules

### ✨ Tidy

- Generate a ‘tidier’ form of the `multiqc_data.json` MultiQC summary
  file. TSV and Parquet outputs are generated. For the TSV file, each
  row corresponds to a single sample, and each column corresponds to a
  single quality metric/variable. See the [CLI](#cli) section below for
  options.

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
    dracarys.R 0.4.0

    #------- Tidy -------#


    $ dracarys.R tidy --help
    usage: dracarys.R tidy [-h] -j JSON -o OUTDIR -p PREFIX [-q]
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
