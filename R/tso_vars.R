tso_snv <- function(snvs) {
  # snvs is an array of snv elements
  snv_info <- function(snv) {
    main <- tibble::tibble(
      chrom = snv[["vcfChromosome"]],
      pos = snv[["vcfPosition"]],
      ref = snv[["vcfRefAllele"]],
      alt = snv[["vcfAltAllele"]],
      af = snv[["vcfVariantFrequency"]],
      qual = snv[["quality"]],
      dp_tot = snv[["totalDepth"]],
      dp_alt = snv[["altAlleleDepth"]]
    )
    # each snv has a single-element array nirvana
    nirv <- snv[["nirvana"]][[1]]
    # nirvana has a transcript array
    txs <- nirv[["transcripts"]]
    get_tx_info <- function(tx) {
      cons <- unlist(tx$consequence) |> paste(collapse = ",")
      tx$consequence <- NULL
      tibble::as_tibble_row(tx) |>
        dplyr::mutate(consequence = cons)
    }
    if (length(txs) > 0) {
      tx_info <- txs |>
        purrr::map_dfr(get_tx_info) |>
        dplyr::bind_rows()
    } else {
      tx_info <- NULL
    }
    dplyr::bind_cols(main, tx_info)
  }
  purrr::map_dfr(snvs, snv_info)
}

tso_cnv <- function(cnvs) {
  purrr::map_dfr(cnvs, tibble::as_tibble)
}
