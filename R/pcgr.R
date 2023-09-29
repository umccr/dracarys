#' PcgrJson R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `pcgr.json.gz` file output from PCGR.
#'
#' @examples
#' \dontrun{
#' x <- "/path/to/pcgr.json.gz"
#' d <- PcgrJsonFile$new(x)
#' d_parsed <- d$read() # or read(d)
#' d$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "both")
#' }
#' @export
PcgrJsonFile <- R6::R6Class(
  "PcgrJsonFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `pcgr.json.gz` file output from PCGR.
    #'
    #' @return List of tibbles.
    read = function() {
      x <- self$path
      j <- read_jsongz_jsonlite(x)
      # l2tib <- function(el) {
      #   purrr::flatten(el) |>
      #     dplyr::bind_rows() |>
      #     dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.character(.)))
      # }
      # dbrel <- j[["metadata"]][["pcgr_db_release"]] |>
      #   purrr::map(l2tib) |>
      #   dplyr::bind_rows(.id = "name_tidy") |>
      #   dplyr::select("name", "name_tidy", "version", "url", "resource_type")
      # handle nulls and rename - see umccr/dracarys#99
      tmb <-
        j[["content"]][["tmb"]][["variant_statistic"]] %||%
        j[["content"]][["tmb"]][["v_stat"]] %||%
        list(tmb_estimate = NA, n_tmb = NA)
      tmb <- purrr::flatten(tmb) |>
        tibble::as_tibble_row() |>
        dplyr::select("tmb_estimate", "n_tmb")
      msi <- j[["content"]][["msi"]][["prediction"]][["msi_stats"]]
      # handle nulls
      msi <- msi %||% list(fracIndels = NA, predicted_class = NA)
      msi <- purrr::flatten(msi) |>
        tibble::as_tibble_row() |>
        dplyr::select("fracIndels", "predicted_class")
      metrics <- dplyr::bind_cols(msi, tmb)
      list(
        # using list in case we want other data as well
        metrics = metrics
      )
    },

    #' @description
    #' Writes a tidy version of the `pcgr.json.gz` file output from PCGR.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s) (one of 'tsv' (def.),
    #' 'parquet', 'both').
    write = function(d, out_dir, prefix, out_format = "tsv") {
      prefix <- file.path(out_dir, prefix)
      p <- glue("{prefix}_pcgr")
      l <- list(
        meta = list(
          obj = d[["metrics"]],
          pref = glue("{p}_metrics")
        )
      )
      purrr::map(l, function(k) {
        write_dracarys(obj = k[["obj"]], prefix = k[["pref"]], out_format = out_format)
      })
    }
  )
)

#' PcgrTiersFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `pcgr.snvs_indels.tiers.tsv` file output from PCGR.
#'
#' @examples
#' \dontrun{
#' x <- "/path/to/pcgr.snvs_indels.tiers.tsv"
#' d <- PcgrTiersFile$new(x)
#' d_parsed <- d$read() # or read(d)
#' d$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "both")
#' }
#' @export
PcgrTiersFile <- R6::R6Class(
  "PcgrTiersFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `pcgr.snvs_indels.tiers.tsv` file output from PCGR.
    #'
    #' @return List of tibbles.
    read = function() {
      x <- self$path
      ct <- readr::cols(
        CHROM = "c", POS = "i", REF = "c", ALT = "c", GENOMIC_CHANGE = "c",
        GENOME_VERSION = "c", VCF_SAMPLE_ID = "c", VARIANT_CLASS = "c",
        SYMBOL = "c", GENE_NAME = "c", CCDS = "c", CANONICAL = "c",
        ENTREZ_ID = "d", UNIPROT_ID = "c", ENSEMBL_TRANSCRIPT_ID = "c",
        ENSEMBL_GENE_ID = "c", REFSEQ_MRNA = "c", ONCOSCORE = "d",
        ONCOGENE = "l", TUMOR_SUPPRESSOR = "l", ONCOGENE_EVIDENCE = "c",
        TUMOR_SUPPRESSOR_EVIDENCE = "c", DISGENET_CUI = "c",
        DISGENET_TERMS = "c", CONSEQUENCE = "c", PROTEIN_CHANGE = "c",
        PROTEIN_DOMAIN = "c", CODING_STATUS = "c", EXONIC_STATUS = "c",
        CDS_CHANGE = "c", HGVSp = "c", HGVSc = "c", EFFECT_PREDICTIONS = "c",
        MUTATION_HOTSPOT = "c", MUTATION_HOTSPOT_TRANSCRIPT = "c",
        MUTATION_HOTSPOT_CANCERTYPE = "c", PUTATIVE_DRIVER_MUTATION = "l",
        CHASMPLUS_DRIVER = "c", CHASMPLUS_TTYPE = "c", VEP_ALL_CSQ = "c",
        DBSNPRSID = "c", COSMIC_MUTATION_ID = "c", TCGA_PANCANCER_COUNT = "d",
        TCGA_FREQUENCY = "c", ICGC_PCAWG_OCCURRENCE = "c",
        CHEMBL_COMPOUND_ID = "c", CHEMBL_COMPOUND_TERMS = "c",
        SIMPLEREPEATS_HIT = "l", WINMASKER_HIT = "l", OPENTARGETS_RANK = "d",
        CLINVAR = "c", CLINVAR_CLNSIG = "c", GLOBAL_AF_GNOMAD = "d",
        GLOBAL_AF_1KG = "d", CALL_CONFIDENCE = "l", DP_TUMOR = "d",
        AF_TUMOR = "d", DP_CONTROL = "l", AF_CONTROL = "l", TIER = "c",
        TIER_DESCRIPTION = "c"
      )
      readr::read_tsv(x, col_types = ct)
    },

    #' @description
    #' Writes a tidy version of the `pcgr.snvs_indels.tiers.tsv` file output from PCGR.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s) (one of 'tsv' (def.),
    #' 'parquet', 'both').
    write = function(d, out_dir, prefix, out_format = "tsv") {
      prefix <- file.path(out_dir, prefix)
      prefix2 <- glue("{prefix}_tiers")
      write_dracarys(obj = d, prefix = prefix2, out_format = out_format)
    }
  )
)
