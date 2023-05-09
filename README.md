
- [🔥 dracarys - UMCCR Workflow
  Tidying](#-dracarys---umccr-workflow-tidying)
  - [🏆 Aim](#-aim)
  - [🍕 Installation](#-installation)
  - [✨ Supported Workflows](#-supported-workflows)
  - [🌀 CLI](#-cli)
  - [🚕 Running](#-running)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# 🔥 dracarys - UMCCR Workflow Tidying

![](https://emojis.slackmojis.com/emojis/images/1643515659/16823/flying_dragon.gif?1643515659 "Dragon Flying")

- Docs: <https://umccr.github.io/dracarys/>

[![Conda
install](https://anaconda.org/umccr/r-dracarys/badges/version.svg)](https://anaconda.org/umccr/r-dracarys)
[![Conda
install](https://anaconda.org/umccr/r-dracarys/badges/latest_release_date.svg)](https://anaconda.org/umccr/r-dracarys)

## 🏆 Aim

Given a [ICA
GDS](https://developer.illumina.com/illumina-connected-analytics) or
local directory with results from a UMCCR workflow, {dracarys} will grab
files of interest and transform them into ‘tidier’ structures for output
into TSV/Parquet format for downstream ingestion into a database/data
lake. See supported [workflows](#supported-workflows),
[running](#running) examples, and [CLI](#cli) options in the sections
below.

## 🍕 Installation

<details>
<summary>
R
</summary>

``` r
remotes::install_github("umccr/dracarys@vX.X.X") # for vX.X.X Release/Tag
```

</details>
<details>
<summary>
Conda
</summary>

- Linux & MacOS (non-M1)

``` bash
mamba create \
  -n dracarys_env \
  -c umccr -c bioconda -c conda-forge \
  r-dracarys==X.X.X

conda activate dracarys_env
```

- MacOS M1

``` bash
CONDA_SUBDIR=osx-64 \
  mamba create \
  -n dracarys_env \
  -c umccr -c bioconda -c conda-forge \
  r-dracarys==X.X.X

conda activate dracarys_env
```

</details>
<details>
<summary>
Docker
</summary>

``` bash
docker pull --platform linux/amd64 ghcr.io/umccr/dracarys:X.X.X
```

</details>

<a name="supported-workflows"></a>

## ✨ Supported Workflows

| Workflow   | Description                                                                                                                                          |
|------------|------------------------------------------------------------------------------------------------------------------------------------------------------|
| 🐉 DRAGEN  | UMCCR [DRAGEN DNA/RNA](https://support-docs.illumina.com/SW/DRAGEN_v40/Content/SW/DRAGEN/GPipelineIntro_fDG.htm) workflows                           |
| 🌶 cttso    | UMCCR [ctDNA TSO500](https://support-docs.illumina.com/SW/DRAGEN_TSO500_ctDNA_v2.1/Content/SW/TSO500/WorkflowDiagram_appT500ctDNAlocal.htm) workflow |
| 🌈 MultiQC | [MultiQC](https://multiqc.info/) QC workflow summaries                                                                                               |

See which output files from these workflows are supported in [Supported
Files](https://umccr.github.io/dracarys/articles/files.html).

<a name="cli"></a>

## 🌀 CLI

A `dracarys.R` command line interface is available for convenience.

- If you’re using the conda package, the `dracarys.R` command will
  already be available inside the activated conda environment.
- If you’re *not* using the conda package, you need to export the
  `dracarys/inst/cli/` directory to your `PATH` in order to use
  `dracarys.R`.

``` bash
dracarys_cli=$(Rscript -e 'x = system.file("cli", package = "dracarys"); cat(x, "\n")' | xargs)
export PATH="${dracarys_cli}:${PATH}"
```

    dracarys.R --version
    dracarys.R 0.8.0

    #-----------------------------------#
    dracarys.R --help
    usage: dracarys.R [-h] [-v] {tidy} ...

    🐉 DRAGEN Output Post-Processing 🔥

    positional arguments:
      {tidy}         sub-command help
        tidy         Tidy UMCCR Workflow Outputs

    options:
      -h, --help     show this help message and exit
      -v, --version  show program's version number and exit

    #-----------------------------------#
    #------- Tidy ----------------------#
    dracarys.R tidy --help
    usage: dracarys.R tidy [-h] -i IN_DIR -o OUT_DIR -p PREFIX [-t TOKEN]
                           [-g GDS_LOCAL_DIR] [-f {tsv,parquet,both}] [-n] [-q]

    options:
      -h, --help            show this help message and exit
      -i IN_DIR, --in_dir IN_DIR
                            ⛄️ Directory with untidy UMCCR workflow results. Can
                            be GDS or local.
      -o OUT_DIR, --out_dir OUT_DIR
                            🔥 Directory to output tidy results.
      -p PREFIX, --prefix PREFIX
                            🎻 Prefix string used for all results.
      -t TOKEN, --token TOKEN
                            🙈 ICA access token. Default: ICA_ACCESS_TOKEN env var.
      -g GDS_LOCAL_DIR, --gds_local_dir GDS_LOCAL_DIR
                            📥 If input is a GDS directory, download the
                            recognisable files to this directory. Default:
                            '<out_dir>/dracarys_gds_sync'.
      -f {tsv,parquet,both}, --format {tsv,parquet,both}
                            🎨 Format of output. Default: tsv.
      -n, --dryrun          🐫 Dry run - just show files to be tidied.
      -q, --quiet           😴 Shush all the logs.

<a name="running"></a>

## 🚕 Running

{dracarys} takes as input (`--in_dir`) a directory with results from one
of the UMCCR [workflows](#supported-workflows). It will recursively scan
that directory for [supported
files](https://umccr.github.io/dracarys/articles/files.html), download
those into a local directory (`--gds_local_dir`), and then it will
parse, transform and write the tidied versions into the specified output
directory (`--out_dir`). A prefix (`--prefix`) is prepended to each of
the tidied files. The output file format (`--format`) can be tsv,
parquet, or both. To get just a list of supported files within the
specified input directory, use the `-n (--dryrun)` option.

<details>
<summary>
R
</summary>

``` r
# help(umccr_tidy)
in_dir <- "gds://path/to/subjectX_multiqc_data/"
out_dir <- tempdir()
prefix <- "subjectX"
umccr_tidy(in_dir = in_dir, out_dir = out_dir, prefix = prefix)
```

</details>
<details>
<summary>
Mac/Linux
</summary>

From within an activated conda environment or a shell with the
`dracarys.R` CLI available:

``` bash
dracarys.R tidy \
      -i gds://path/to/subjectX_multiqc_data/ \
      -o local_output_dir \
      -p subjectX_prefix
```

</details>
<details>
<summary>
Docker
</summary>

``` bash
docker container run \
  -v $(PWD):/mount1 \
  --platform=linux/amd64 \
  --env "ICA_ACCESS_TOKEN" \
  --rm -it \
  ghcr.io/umccr/dracarys:X.X.X \
    dracarys.R tidy \
      -i gds://path/to/subjectX_multiqc_data/ \
      -o /mount1/output_dir \
      -p subjectX_prefix
```

</details>
