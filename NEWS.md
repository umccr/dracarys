# dracarys 0.11.0 (2023-08-22)

[0.10.0 - 0.11.0 diff](https://github.com/umccr/dracarys/compare/v0.10.0...v0.11.0)

# dracarys 0.10.0 (2023-06-22)

[0.9.0 - 0.10.0 diff](https://github.com/umccr/dracarys/compare/v0.9.0...v0.10.0)

# dracarys 0.9.0 (2023-05-10)

[0.8.0 - 0.9.0 diff](https://github.com/umccr/dracarys/compare/v0.8.0...v0.9.0)

# dracarys 0.8.0 (2023-02-24)

[0.7.0 - 0.8.0 diff](https://github.com/umccr/dracarys/compare/v0.7.0...v0.8.0)

# dracarys 0.7.0 (2023-02-15)

[0.6.0 - 0.7.0 diff](https://github.com/umccr/dracarys/compare/v0.6.0...v0.7.0)

Mostly added use cases based on requests from the curation team related to
umccrise (CHORD, HRDetect, QC summary), PCGR, and MultiQC.

- :star: Added :
  - R functions:
    - `umccr_tidy`: workhorse that handles the tidying and method dispatch.
    - `gds_file_presignedurl`
  - R6 classes:
    - `PcgrJsonFile`
    - `PcgrTiersFile`
    - `UmChordTsvFile`
    - `UmHrdetectTsvFile`
    - `UmQcSumFile`
    - `UmSigsSnvFile`
    - `MultiqcFile`
- :wrench:
  - Updated README with unified `dracarys tidy` CLI.
  - Added a `fun` column to `FILE_REGEX` for handling method dispatch.
  - GH Actions: Remove linux/arm64 from docker since conda pkgs don't generally have ARM64 equivalents.
  - GH Actions: use "miniforge-variant: Mambaforge"
    (see [this issue](https://github.com/conda-incubator/setup-miniconda/issues/274)).

# dracarys 0.6.0 (2023-01-09)

[0.5.0 - 0.6.0 diff](https://github.com/umccr/dracarys/compare/v0.5.0...v0.6.0)

- :wrench: GH Actions:
  - replace `::set-output`.
  - add linux/arm64 support.
- :wrench: R core:
  - replace select `.data` with quotes.
  - bump Roxygen `7.2.1 -> 7.2.2`
  - add {fs}, {httr}, {jose}, {sessioninfo} dependencies.
  - import `%||%` from {rlang} to specify fallback values in NULL cases.
- :sparkles: pkgdown: change theme to simplex.
- :whale: Docker:
  - bump mambaforge: `4.12.0-2 -> 22.9.0-2`
  - bump conda-lock: `1.0.5 -> 1.3.0`
- :computer: CLI:
  - add support for tidying outputs from the UMCCR TSO500 ctDNA workflow, and
    for generating a HTML file with tidy results from the same workflow.
- :star: Added R functions:
  - `tso_tidy` for tidying outputs from the UMCCR TSO500 ctDNA workflow.
  - `tso_rmd` for generating a HTML file with tidy results from the above workflow.
  - `dr_gds_download` for downloading dracarys-related files from GDS.
  - `gds_files_list` for listing files on GDS via the API.
  - `gds_file_download` and `gds_file_download_api` for downloading files
    from GDS via the ica binary and API (via presigned URLs), respectively.
  - `ica_token_validate` for validating the ICA access token is valid and has not expired.
- :x: Removed R functions:
  - `TsoCombinedVariantOutputFile`: these files were deemed to have less info compared to other output files.
  - `dracarys_tidy_multiqc`: renamed to `dracarys_multiqc`.
  - `mkdir`: use `fs::create_dir` instead.
  - Remove {tibble}, {readr} and {ggplot} multi-imports for R6, since we can get away with just importing a single
    function in one of the R6 classes.

# dracarys 0.5.0 (2022-09-28)

[0.4.0 - 0.5.0 diff](https://github.com/umccr/dracarys/compare/v0.4.0...v0.5.0)

- MultiQC: update column mappings ([pr15](https://github.com/umccr/dracarys/pull/15), [pr16](https://github.com/umccr/dracarys/pull/16)).
  - move map to separate TSV
- CLI: add option for output format (tsv, parquet, or both) ([pr18](https://github.com/umccr/dracarys/pull/18)).
- new contributors: [@victorskl](https://github.com/victorskl)

# dracarys 0.4.0 (2022-09-12)

[0.3.0 - 0.4.0 diff](https://github.com/umccr/dracarys/compare/v0.3.0...v0.4.0)

- :star: support for DRAGEN TSO500 ctdna output ([pr14](https://github.com/umccr/dracarys/pull/14)).
  - also add Quarto HTML report template

# dracarys 0.3.0 (2022-08-28)

[0.2.0 - 0.3.0 diff](https://github.com/umccr/dracarys/compare/v0.2.0...v0.3.0)

- :star: MultiQC: support for DRAGEN ctdna output ([pr13](https://github.com/umccr/dracarys/pull/13)).

# dracarys 0.2.0 (2022-07-23)

[0.1.0 - 0.2.0 diff](https://github.com/umccr/dracarys/compare/v0.1.0...v0.2.0)

- Add `MULTIQC_COLUMNS` tibble that maps the raw metric name to a cleaner name
  ([pr12](https://github.com/umccr/dracarys/pull/12)).
- Support for more MultiQC JSONs (from bcbio-wts and bcbio-wgs).

# dracarys 0.1.0 (2022-07-11)

- Initial release of dracarys.
- Support for [MultiQC](https://github.com/ewels/MultiQC) JSON tidying.
- Add conda, conda-lock, Docker support ([pr11](https://github.com/umccr/dracarys/pull/11)).
- CLI support for `dracarys.R tidy`.
  - Initially supporting just MultiQC json input (from dragen-alignment,
    dragen-transcriptome, dragen-tumor-normal, dragen-umccrise,
    and bcbio-umccrise).
