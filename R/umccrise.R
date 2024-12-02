#' Wf_umccrise R6 Class
#'
#' @description
#' Reads and writes tidy versions of files from the `umccrise` workflow
#'
#' @examples
#' \dontrun{
#'
#' #---- LOCAL ----#
#' SubjectID <- "SBJ04662"
#' SampleID_tumor <- "PRJ240647"
#' SampleID_normal <- "PRJ240646"
#' p1_local <- "~/icav1/g/production/analysis_data"
#' p <- file.path(normalizePath(p1_local), "SBJ04662/umccrise/20240302e66750fe/L2400240__L2400239")
#' um1 <- Wf_umccrise$new(
#'   path = p, SubjectID = SubjectID,
#'   SampleID_tumor = SampleID_tumor, SampleID_normal = SampleID_normal
#' )
#' um1$list_files(max_files = 10)
#' um1$list_files_filter_relevant()
#' d <- um1$download_files(max_files = 1000, dryrun = F)
#' d_tidy <- um1$tidy_files(d)
#' d_write <- um1$write(
#'   d_tidy,
#'   outdir = file.path(p, "dracarys_tidy"),
#'   prefix = glue("{SubjectID}__{SampleID_tumor}"),
#'   format = "tsv"
#' )
#'
#' #---- GDS ----#
#' SubjectID <- "SBJ04662"
#' SampleID_tumor <- "PRJ240647"
#' SampleID_normal <- "PRJ240646"
#' p1_gds <- "gds://production/analysis_data"
#' p <- file.path(p1_gds, "SBJ04662/umccrise/20240302e66750fe/L2400240__L2400239")
#' outdir <- file.path(sub("gds:/", "~/icav1/g", p))
#' token <- Sys.getenv("ICA_ACCESS_TOKEN")
#' um2 <- Wf_umccrise$new(
#'   path = p, SubjectID = SubjectID,
#'   SampleID_tumor = SampleID_tumor, SampleID_normal = SampleID_normal
#' )
#' um2$list_files(max_files = 8)
#' um2$list_files_filter_relevant(ica_token = token, max_files = 1000)
#' d <- um2$download_files(
#'   outdir = outdir, ica_token = token,
#'   max_files = 1000, dryrun = F
#' )
#' d_tidy <- um2$tidy_files(d)
#' d_write <- um2$write(
#'   d_tidy,
#'   outdir = file.path(outdir, "dracarys_tidy"),
#'   prefix = glue("{SubjectID}__{SampleID_tumor}"),
#'   format = "tsv"
#' )
#' }
#'
#' @export
Wf_umccrise <- R6::R6Class(
  "Wf_umccrise",
  inherit = Wf,
  public = list(
    #' @field SubjectID The SubjectID of the sample.
    #' @field SampleID_tumor The SampleID of the tumor sample.
    #' @field SampleID_normal The SampleID of the normal sample.
    SubjectID = NULL,
    SampleID_tumor = NULL,
    SampleID_normal = NULL,
    #' @description Create a new Wf_umccrise object.
    #' @param path Path to directory with raw workflow results (from GDS, S3, or
    #' local filesystem).
    #' @param SubjectID The SubjectID of the sample.
    #' @param SampleID_tumor The SampleID of the tumor sample.
    #' @param SampleID_normal The SampleID of the normal sample.
    initialize = function(path = NULL, SubjectID = NULL,
                          SampleID_tumor = NULL, SampleID_normal = NULL) {
      wname <- "umccrise"
      pref <- glue("{SubjectID}__{SampleID_tumor}")
      crep <- "cancer_report_tables"
      regexes <- tibble::tribble(
        ~regex, ~fun,
        glue("{path}/{pref}/{crep}/hrd/{pref}-chord\\.tsv\\.gz$"), "read_hrdChord",
        glue("{path}/{pref}/{crep}/hrd/{pref}-hrdetect\\.tsv\\.gz$"), "read_hrdHrdetect",
        glue("{path}/{pref}/{crep}/sigs/{pref}-snv_2015\\.tsv\\.gz$"), "read_sigsTsv",
        glue("{path}/{pref}/{crep}/sigs/{pref}-snv_2020\\.tsv\\.gz$"), "read_sigsTsv",
        glue("{path}/{pref}/{crep}/sigs/{pref}-dbs\\.tsv\\.gz$"), "read_sigsTsv",
        glue("{path}/{pref}/{crep}/sigs/{pref}-indel\\.tsv\\.gz$"), "read_sigsTsv",
        glue("{path}/{pref}/{crep}/{pref}-qc_summary\\.tsv\\.gz$"), "read_qcSum",
        glue("{path}/{pref}/{pref}-multiqc_report_data/multiqc_conpair\\.txt$"), "read_conpair",
        glue("{path}/{pref}/purple/{pref}\\.purple\\.cnv\\.gene\\.tsv$"), "DOWNLOAD_ONLY-purplegene",
        glue("{path}/{pref}/small_variants/{pref}-somatic-PASS\\.vcf\\.gz$"), "DOWNLOAD_ONLY-smlvfiltvcf",
        glue("{path}/{pref}/small_variants/{pref}-somatic-PASS\\.vcf\\.gz\\.tbi$"), "DOWNLOAD_ONLY-smlvfiltvcfi",
        glue("{path}/work/{pref}/pcgr/{pref}-somatic\\.pcgr\\.json\\.gz$"), "read_pcgrJson",
        glue("{path}/work/{pref}/pcgr/{pref}-somatic\\.pcgr\\.pass\\.vcf\\.gz$"), "DOWNLOAD_ONLY-pcgrvcf",
        glue("{path}/work/{pref}/pcgr/{pref}-somatic\\.pcgr\\.pass\\.vcf\\.gz\\.tbi$"), "DOWNLOAD_ONLY-pcgrvcfi",
        glue("{path}/work/{pref}/pcgr/{pref}-somatic\\.pcgr\\.snvs_indels\\.tiers\\.tsv$"), "DOWNLOAD_ONLY-pcgrtiers",
        glue("{path}/work/{pref}/cpsr/{pref}-normal\\.cpsr\\.vcf\\.gz$"), "DOWNLOAD_ONLY-cpsrvcf",
        glue("{path}/work/{pref}/cpsr/{pref}-normal\\.cpsr\\.vcf\\.gz\\.tbi$"), "DOWNLOAD_ONLY-cpsrvcfi"
      )
      super$initialize(path = path, wname = wname, regexes = regexes)
      self$SubjectID <- SubjectID
      self$SampleID_tumor <- SampleID_tumor
      self$SampleID_normal <- SampleID_normal
    },
    #' @description Print details about the Workflow.
    #' @param ... (ignored).
    print = function(...) {
      res <- tibble::tribble(
        ~var, ~value,
        "path", private$.path,
        "wname", private$.wname,
        "filesystem", private$.filesystem,
        "nregexes", as.character(nrow(private$.regexes)),
        "SubjectID", self$SubjectID,
        "SampleID_tumor", self$SampleID_tumor,
        "SampleID_normal", self$SampleID_normal
      )
      print(res)
      invisible(self)
    },
    #' @description Read `pcgr.json.gz` file.
    #' @param x Path to file.
    read_pcgrJson = function(x) {
      dat <- pcgr_json_read(x)
      tibble::tibble(name = "pcgrjson", data = list(dat))
    },
    #' @description Read `chord.tsv.gz` cancer report file.
    #' @param x Path to file.
    read_hrdChord = function(x) {
      ct <- readr::cols_only(
        p_hrd = "d",
        hr_status = "c",
        hrd_type = "c",
        p_BRCA1 = "d",
        p_BRCA2 = "d"
      )
      dat <- read_tsvgz(x, col_types = ct)
      tibble::tibble(name = "hrdchord", data = list(dat[]))
    },
    #' @description Read `hrdetect.tsv.gz` cancer report file.
    #' @param x Path to file.
    read_hrdHrdetect = function(x) {
      ct <- readr::cols(
        .default = "d",
        sample = "c"
      )
      dat <- read_tsvgz(x, col_types = ct) |>
        dplyr::select(-c("sample"))
      tibble::tibble(name = "hrdhrdetect", data = list(dat[]))
    },
    #' @description Read signature cancer report file.
    #' @param x Path to file.
    read_sigsTsv = function(x) {
      .sigsSuffix <- function(x) {
        x <- basename(x)
        dplyr::case_when(
          grepl("-dbs", x) ~ "dbs",
          grepl("-indel", x) ~ "ind",
          grepl("-snv_2015", x) ~ "snv2015",
          grepl("-snv_2020", x) ~ "snv2020",
          .default = ""
        )
      }
      suffix <- .sigsSuffix(x)
      ct <- readr::cols(
        .default = "d",
        Signature = "c"
      )
      dat <- read_tsvgz(x, col_types = ct)
      tibble::tibble(name = glue("sigs_{suffix}"), data = list(dat[]))
    },
    #' @description Read `qc_summary.tsv.gz` cancer report file.
    #' @param x Path to file.
    read_qcSum = function(x) {
      d <- read_tsvgz(x, col_types = readr::cols(.default = "c"))
      dat <- d |>
        dplyr::select("variable", "value") |>
        tidyr::pivot_wider(names_from = "variable", values_from = "value") |>
        dplyr::rename(MSI_mb_tmp = "MSI (indels/Mb)") |>
        dplyr::mutate(
          purity_hmf = sub("(.*) \\(.*\\)", "\\1", .data$Purity) |> as.numeric(),
          ploidy_hmf = sub("(.*) \\(.*\\)", "\\1", .data$Ploidy) |> as.numeric(),
          hrd_chord = sub("CHORD: (.*); HRDetect: (.*)", "\\1", .data$HRD) |> as.numeric(),
          hrd_hrdetect = sub("CHORD: (.*); HRDetect: (.*)", "\\2", .data$HRD),
          # handle HRDetect NA
          hrd_hrdetect = ifelse(.data$hrd_hrdetect == "NA", NA_real_, as.numeric(.data$hrd_hrdetect)),
          msi_mb_hmf = sub(".* \\((.*)\\)", "\\1", .data$MSI_mb_tmp) |> as.numeric(),
          contamination_hmf = as.numeric(.data$Contamination),
          deleted_genes_hmf = as.numeric(.data$DeletedGenes),
          msi_hmf = sub("(.*) \\(.*\\)", "\\1", .data$MSI_mb_tmp),
          tmb_hmf = sub("(.*) \\(.*\\)", "\\1", .data$TMB) |> as.numeric(),
          tml_hmf = sub("(.*) \\(.*\\)", "\\1", .data$TML) |> as.numeric(),
          hypermutated = ifelse("Hypermutated" %in% d$variable, .data[["Hypermutated"]], NA) |> as.character(),
          bpi_enabled = ifelse("BPI Enabled" %in% d$variable, .data[["BPI Enabled"]], NA) |> as.character(),
        ) |>
        dplyr::select(
          qc_status_hmf = "QC_Status",
          sex_hmf = "Gender",
          "purity_hmf", "ploidy_hmf", "msi_hmf", "msi_mb_hmf",
          "hrd_chord", "hrd_hrdetect", "contamination_hmf",
          "deleted_genes_hmf", "tmb_hmf", "tml_hmf",
          wgd_hmf = "WGD",
          "hypermutated", "bpi_enabled"
        )
      tibble::tibble(name = glue("qcsum"), data = list(dat[]))
    },
    #' @description Read multiqc_conpair.txt file.
    #' @param x Path to file.
    read_conpair = function(x) {
      um_ref_samples <- c("Alice", "Bob", "Chen", "Elon", "Dakota")
      um_ref_samples <- paste0(um_ref_samples, rep(c("_T", "_B", ""), each = length(um_ref_samples)))
      cnames <- list(
        old = c(
          "Sample", "concordance_concordance", "concordance_used_markers",
          "concordance_total_markers", "concordance_marker_threshold",
          "concordance_min_mapping_quality", "concordance_min_base_quality",
          "contamination"
        ),
        new = c(
          "sampleid", "contamination", "concordance", "markers_used",
          "markers_total", "marker_threshold",
          "mapq_min", "baseq_min"
        )
      )
      ctypes <- list(
        old = c("cddddddd"),
        new = c("cddddddd")
      )
      if (!file.exists(x)) {
        return(empty_tbl(cnames$new, ctypes$new))
      }
      d1 <- readr::read_tsv(x, col_types = readr::cols(.default = "d", Sample = "c"))
      assertthat::assert_that(all(colnames(d1) == cnames$old))
      dat <- d1 |>
        dplyr::filter(!.data$Sample %in% um_ref_samples) |>
        dplyr::relocate("contamination", .after = "Sample") |>
        rlang::set_names(cnames$new)
      tibble::tibble(name = glue("conpair"), data = list(dat[]))
    }
  ) # end public
)

#' umccrise Download Tidy and Write
#'
#' Downloads files from the `umccrise` workflow and writes them in a tidy format.
#'
#' @param path Path to directory with raw workflow results (from GDS, S3, or
#' local filesystem).
#' @param SubjectID The SubjectID of the sample.
#' @param SampleID_tumor The SampleID of the tumor sample.
#' @param SampleID_normal The SampleID of the normal sample.
#' @param outdir Path to output directory.
#' @param format Format of output files.
#' @param max_files Max number of files to list.
#' @param ica_token ICA access token (def: $ICA_ACCESS_TOKEN env var).
#' @param dryrun If TRUE, just list the files that will be downloaded (don't
#' download them).
#' @return Tibble of tidy data as list-cols.
#'
#' @examples
#' \dontrun{
#' SubjectID <- "SBJ04662"
#' SampleID_tumor <- "PRJ240647"
#' SampleID_normal <- "PRJ240646"
#' p1_gds <- "gds://production/analysis_data"
#' p <- file.path(p1_gds, "SBJ04662/umccrise/20240302e66750fe/L2400240__L2400239")
#' outdir <- file.path(sub("gds:/", "~/icav1/g", p))
#' token <- Sys.getenv("ICA_ACCESS_TOKEN")
#' d <- Wf_umccrise_download_tidy_write(
#'   path = p, SubjectID = SubjectID,
#'   SampleID_tumor = SampleID_tumor, SampleID_normal = SampleID_normal,
#'   outdir = outdir,
#'   dryrun = F
#' )
#' }
#' @export
Wf_umccrise_download_tidy_write <- function(path, SubjectID,
                                            SampleID_tumor, SampleID_normal,
                                            outdir, format = "rds", max_files = 1000,
                                            ica_token = Sys.getenv("ICA_ACCESS_TOKEN"),
                                            dryrun = FALSE) {
  um <- Wf_umccrise$new(
    path = path, SubjectID = SubjectID,
    SampleID_tumor = SampleID_tumor, SampleID_normal = SampleID_normal
  )
  d_dl <- um$download_files(
    outdir = outdir, ica_token = ica_token,
    max_files = max_files, dryrun = dryrun
  )
  if (!dryrun) {
    d_tidy <- um$tidy_files(d_dl)
    d_write <- um$write(
      d_tidy,
      outdir = file.path(outdir, "dracarys_tidy"),
      prefix = glue("{SubjectID}__{SampleID_tumor}"),
      format = format
    )
    return(d_write)
  }
  return(d_dl)
}
