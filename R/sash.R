#' Wf_sash R6 Class
#'
#' @description
#' Reads and writes tidy versions of files from the `sash` workflow
#'
#' @examples
#' \dontrun{
#'
#' #---- Local ----#
#' p1 <- "~/s3/org.umccr.data.oncoanalyser/analysis_data/SBJ02862/sash"
#' p2 <- "20240830ece6b0b7/L2201449_L2201450"
#' p <- normalizePath(file.path(p1, p2))
#' SubjectID <- "SBJ02862"
#' SampleID_tumor <- "PRJ222112"
#' SampleID_normal <- "PRJ222114"
#' prefix <- glue("{SubjectID}_{SampleID_tumor}")
#' s1 <- Wf_sash$new(
#'   path = p, SubjectID = SubjectID,
#'   SampleID_tumor = SampleID_tumor, SampleID_normal = SampleID_normal
#' )
#' #-- test regexes active binding
#' counts1 <- glue(
#'   "{p}/{prefix}/smlv_somatic/report/",
#'   "{SampleID_tumor}\\.somatic\\.variant_counts_process\\.json$"
#' )
#' regexes1 <- tibble::tribble(
#'   ~regex, ~fun,
#'   counts1, "read_smlvSomCounts"
#' )
#' s1$regexes <- regexes1
#' s1$list_files(max_files = 20)
#' s1$list_files_filter_relevant(max_files = 300)
#' d <- s1$download_files(max_files = 1000, dryrun = F)
#' d_tidy <- s1$tidy_files(d)
#' d_write <- s1$write(
#'   d_tidy,
#'   outdir = file.path(p, "dracarys_tidy"),
#'   prefix = glue("{SubjectID}_{SampleID_tumor}"),
#'   format = "tsv"
#' )
#'
#' #---- S3 ----#
#' p1 <- "s3://org.umccr.data.oncoanalyser/analysis_data/SBJ02862/sash"
#' p2 <- "20240830ece6b0b7/L2201449_L2201450"
#' p <- file.path(p1, p2)
#' SubjectID <- "SBJ02862"
#' SampleID_tumor <- "PRJ222112"
#' SampleID_normal <- "PRJ222114"
#' s1 <- Wf_sash$new(
#'   path = p, SubjectID = SubjectID,
#'   SampleID_tumor = SampleID_tumor, SampleID_normal = SampleID_normal
#' )
#' s1$list_files(max_files = 20)
#' s1$list_files_filter_relevant()
#' outdir <- sub("s3:/", "~/s3", p)
#' d <- s1$download_files(outdir = outdir, max_files = 1000, dryrun = F)
#' d_tidy <- s1$tidy_files(d)
#' d_write <- s1$write(
#'   d_tidy,
#'   outdir = file.path(p, "dracarys_tidy"),
#'   prefix = glue("{SubjectID}__{SampleID_tumor}"),
#'   format = "tsv"
#' )
#' }
#'
#' @export
Wf_sash <- R6::R6Class(
  "Wf_sash",
  inherit = Wf,
  public = list(
    #' @field SubjectID The SubjectID of the sample.
    #' @field SampleID_tumor The SampleID of the tumor sample.
    #' @field SampleID_normal The SampleID of the normal sample.
    SubjectID = NULL,
    SampleID_tumor = NULL,
    SampleID_normal = NULL,
    #' @description Create a new Wf_sash object.
    #' @param path Path to directory with raw workflow results (from S3 or
    #' local filesystem).
    #' @param SubjectID The SubjectID of the sample.
    #' @param SampleID_tumor The SampleID of the tumor sample.
    #' @param SampleID_normal The SampleID of the normal sample.
    initialize = function(
      path = NULL,
      SubjectID = NULL,
      SampleID_tumor = NULL,
      SampleID_normal = NULL
    ) {
      wname <- "sash"
      pref <- glue("{SubjectID}_{SampleID_tumor}")
      crep <- "cancer_report/cancer_report_tables"
      # fmt: skip
      regexes <- tibble::tribble(
        ~regex, ~fun,
        glue("{path}/{pref}/{crep}/hrd/{pref}-chord\\.tsv\\.gz$"), "read_hrdChord",
        glue("{path}/{pref}/{crep}/hrd/{pref}-hrdetect\\.tsv\\.gz$"), "read_hrdHrdetect",
        glue("{path}/{pref}/{crep}/hrd/{pref}-dragen\\.tsv\\.gz$"), "read_hrdDragen",
        glue("{path}/{pref}/{crep}/sigs/{pref}-snv_2015\\.tsv\\.gz$"), "read_sigsTsv",
        glue("{path}/{pref}/{crep}/sigs/{pref}-snv_2020\\.tsv\\.gz$"), "read_sigsTsv",
        glue("{path}/{pref}/{crep}/sigs/{pref}-dbs\\.tsv\\.gz$"), "read_sigsTsv",
        glue("{path}/{pref}/{crep}/sigs/{pref}-indel\\.tsv\\.gz$"), "read_sigsTsv",
        glue("{path}/{pref}/{crep}/{pref}-qc_summary\\.tsv\\.gz$"), "read_qcSum",
        glue("{path}/{pref}/purple/{SampleID_tumor}\\.purple\\.cnv\\.gene\\.tsv$"), "DOWNLOAD_ONLY-purplegene",
        glue("{path}/{pref}/smlv_somatic/filter/{SampleID_tumor}\\.pass\\.vcf\\.gz$"), "DOWNLOAD_ONLY-smlvfiltvcf",
        glue("{path}/{pref}/smlv_somatic/filter/{SampleID_tumor}\\.pass\\.vcf\\.gz\\.tbi$"), "DOWNLOAD_ONLY-smlvfiltvcfi",
        glue("{path}/{pref}/smlv_somatic/report/pcgr/{SampleID_tumor}\\.pcgr_acmg\\.grch38\\.json\\.gz$"), "read_pcgrJson",
        glue("{path}/{pref}/smlv_somatic/report/pcgr/{SampleID_tumor}\\.pcgr_acmg\\.grch38\\.vcf\\.gz$"), "DOWNLOAD_ONLY-pcgrvcf",
        glue("{path}/{pref}/smlv_somatic/report/pcgr/{SampleID_tumor}\\.pcgr_acmg\\.grch38\\.vcf\\.gz\\.tbi$"), "DOWNLOAD_ONLY-pcgrvcfi",
        glue("{path}/{pref}/smlv_somatic/report/pcgr/{SampleID_tumor}\\.pcgr_acmg\\.grch38\\.snvs_indels\\.tiers\\.tsv$"), "DOWNLOAD_ONLY-pcgrtiers",
        # glue("{path}/{pref}/smlv_somatic/report/{SampleID_tumor}\\.somatic\\.variant_counts_process\\.json$"), "smlvSomCounts",
        glue("{path}/{pref}/smlv_germline/report/cpsr/{SampleID_normal}\\.cpsr\\.grch38\\.vcf\\.gz$"), "DOWNLOAD_ONLY-cpsrvcf",
        glue("{path}/{pref}/smlv_germline/report/cpsr/{SampleID_normal}\\.cpsr\\.grch38\\.vcf\\.gz\\.tbi$"), "DOWNLOAD_ONLY-cpsrvcfi",
      )
      super$initialize(path = path, wname = wname, regexes = regexes)
      self$SubjectID <- SubjectID
      self$SampleID_tumor <- SampleID_tumor
      self$SampleID_normal <- SampleID_normal
    },
    #' @description Print details about the Workflow.
    #' @param ... (ignored).
    print = function(...) {
      # fmt: skip
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
    #' @description Read `somatic.variant_counts_process.json` file.
    #' @param x Path to file.
    read_smlvSomCounts = function(x) {
      dat <- jsonlite::read_json(x) |>
        tibble::as_tibble_row()
      tibble::tibble(name = "smlvsomcounts", data = list(dat[]))
    },
    #' @description Read `pcgr.json.gz` file.
    #' @param x Path to file.
    read_pcgrJson = function(x) {
      dat <- pcgr_json_read(x)
      tibble::tibble(name = "pcgrjson", data = list(dat[]))
    },
    #' @description Read `dragen.tsv.gz` cancer report hrd file.
    #' @param x Path to file.
    read_hrdDragen = function(x) {
      ct <- readr::cols(.default = "d", Sample = "c")
      dat <- read_tsvgz(x, col_types = ct)
      tibble::tibble(name = "hrddragen", data = list(dat[]))
    },
    #' @description Read `chord.tsv.gz` cancer report hrd file.
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
    #' @description Read `hrdetect.tsv.gz` cancer report hrd file.
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
          purity_hmf = sub("(.*) \\(.*\\)", "\\1", .data$Purity) |>
            as.numeric(),
          ploidy_hmf = sub("(.*) \\(.*\\)", "\\1", .data$Ploidy) |>
            as.numeric(),
          msi_mb_hmf = sub(".* \\((.*)\\)", "\\1", .data$MSI_mb_tmp) |>
            as.numeric(),
          contamination_hmf = as.numeric(.data$Contamination),
          deleted_genes_hmf = as.numeric(.data$DeletedGenes),
          msi_hmf = sub("(.*) \\(.*\\)", "\\1", .data$MSI_mb_tmp),
          tmb_hmf = sub("(.*) \\(.*\\)", "\\1", .data$TMB) |> as.numeric(),
          tml_hmf = sub("(.*) \\(.*\\)", "\\1", .data$TML) |> as.numeric(),
          hypermutated = ifelse(
            "Hypermutated" %in% d$variable,
            .data[["Hypermutated"]],
            NA
          ) |>
            as.character()
        ) |>
        dplyr::select(
          qc_status_hmf = "QC_Status",
          sex_hmf = "Gender",
          "purity_hmf",
          "ploidy_hmf",
          "msi_hmf",
          "msi_mb_hmf",
          "contamination_hmf",
          "deleted_genes_hmf",
          "tmb_hmf",
          "tml_hmf",
          wgd_hmf = "WGD",
          "hypermutated"
        )
      tibble::tibble(name = glue("qcsum"), data = list(dat[]))
    }
  ) # end public
)

#' sash Download Tidy and Write
#'
#' Downloads files from the `sash` workflow and writes them in a tidy format.
#'
#' @param path Path to directory with raw workflow results (from S3 or
#' local filesystem).
#' @param SubjectID The SubjectID of the sample.
#' @param SampleID_tumor The SampleID of the tumor sample.
#' @param SampleID_normal The SampleID of the normal sample.
#' @param outdir Path to output directory.
#' @param format Format of output files.
#' @param max_files Max number of files to list.
#' @param dryrun If TRUE, just list the files that will be downloaded (don't
#' download them).
#' @param regexes Tibble with file `regex` and `fun`ction to parse it. Use only
#' if you want to override the default regexes used for this workflow.
#'
#'
#' @return List where each element is a tidy tibble of a sash file.
#'
#' @examples
#' \dontrun{
#' SubjectID <- "SBJ03043"
#' SampleID_tumor <- "PRJ230004"
#' p1_gds <- glue("gds://production/analysis_data/{SubjectID}/umccrise")
#' p <- file.path(p1_gds, "20240830ec648f40/L2300064__L2300063")
#' outdir <- file.path(sub("gds:/", "~/icav1/g", p))
#' d <- Wf_sash_download_tidy_write(
#'   path = p, SubjectID = SubjectID, SampleID_tumor = SampleID_tumor,
#'   outdir = outdir,
#'   dryrun = F
#' )
#' }
#' @export
Wf_sash_download_tidy_write <- function(
  path,
  SubjectID,
  SampleID_tumor,
  SampleID_normal,
  outdir,
  format = "rds",
  max_files = 1000,
  regexes = NULL,
  dryrun = FALSE
) {
  s <- Wf_sash$new(
    path = path,
    SubjectID = SubjectID,
    SampleID_tumor = SampleID_tumor,
    SampleID_normal = SampleID_normal
  )
  if (!is.null(regexes)) {
    s$regexes <- regexes
  }
  d_dl <- s$download_files(
    outdir = outdir,
    max_files = max_files,
    dryrun = dryrun
  )
  if (!dryrun) {
    d_tidy <- s$tidy_files(d_dl)
    d_write <- s$write(
      d_tidy,
      outdir = file.path(outdir, "dracarys_tidy"),
      prefix = glue("{SubjectID}__{SampleID_tumor}"),
      format = format
    )
    return(d_write)
  }
  return(d_dl)
}
