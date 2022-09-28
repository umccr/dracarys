# dracarys 0.5.0 (2022-09-28)

- MultiQC: update column mappings ([pr15](https://github.com/umccr/dracarys/pull/15), [pr16](https://github.com/umccr/dracarys/pull/16)).
  - move map to separate TSV
- CLI: add option for output format (tsv, parquet, or both) ([pr17](https://github.com/umccr/dracarys/pull/17)).
- contributors: [@victorskl](https://github.com/victorskl)

# dracarys 0.4.0 (2022-09-12)

- :star: support for DRAGEN TSO500 ctdna output ([pr14](https://github.com/umccr/dracarys/pull/14)).
  - also add Quarto HTML report template

# dracarys 0.3.0 (2022-08-28)

- :star: MultiQC: support for DRAGEN ctdna output ([pr13](https://github.com/umccr/dracarys/pull/13)).

# dracarys 0.2.0 (2022-07-23)

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
