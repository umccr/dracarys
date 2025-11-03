#' Wf_sash R6 Class
#'
#' @description
#' Reads and writes tidy versions of files from the `sash` workflow
#'
#' @examples
#' \dontrun{
#'
#' #---- Local ----#
#' p1 <- "~/s3/project-data-889522050439-ap-southeast-2/byob-icav2"
#' p2 <- "project-wgs-accreditation/analysis/sash/20250919121765d8"
#' p <- normalizePath(file.path(p1, p2))
#' libid_tumor <- "L2300943"
#' libid_normal <- "L2300950"
#' s1 <- Wf_sash$new(path = p, libid_tumor = libid_tumor, libid_normal = libid_normal)
#' s1$list_files(max_files = 20)
#' s1$list_files_filter_relevant(max_files = 300)
#' d <- s1$download_files(max_files = 1000, dryrun = F)
#' d_tidy <- s1$tidy_files(d)
#' d_write <- s1$write(
#'   d_tidy,
#'   outdir = file.path(p, "dracarys_tidy"),
#'   prefix = glue("{libid_normal}_{libid_tumor}"),
#'   format = "tsv"
#' )
#'
#' #---- S3 ----#
#' p1 <- "s3://org.umccr.data.oncoanalyser/analysis_data/SBJ02862/sash"
#' p2 <- "20240830ece6b0b7/L2201449_L2201450"
#' p <- file.path(p1, p2)
#' libid_tumor <- "PRJ222112"
#' libid_normal <- "PRJ222114"
#' s1 <- Wf_sash$new(
#'   path = p,
#'   libid_tumor = libid_tumor, libid_normal = libid_normal
#' )
#' s1$list_files(max_files = 20)
#' s1$list_files_filter_relevant()
#' outdir <- sub("s3:/", "~/s3", p)
#' d <- s1$download_files(outdir = outdir, max_files = 1000, dryrun = F)
#' d_tidy <- s1$tidy_files(d)
#' d_write <- s1$write(
#'   d_tidy,
#'   outdir = file.path(p, "dracarys_tidy"),
#'   prefix = glue("__{libid_tumor}"),
#'   format = "tsv"
#' )
#' }
#'
#' @export
Wf_sash <- R6::R6Class(
  "Wf_sash",
  inherit = Wf,
  public = list(
    #' @field libid_tumor The LibraryID of the tumor sample.
    #' @field libid_normal The LibraryID of the normal sample.
    libid_tumor = NULL,
    libid_normal = NULL,
    #' @description Create a new Wf_sash object.
    #' @param path Path to directory with raw workflow results (from S3 or
    #' local filesystem).
    #' @param libid_tumor The LibraryID of the tumor sample.
    #' @param libid_normal The LibraryID of the normal sample.
    initialize = function(path = NULL, libid_tumor = NULL, libid_normal = NULL) {
      wname <- "sash"
      batch <- glue("{libid_tumor}__{libid_normal}")
      pref <- glue("{batch}_{libid_tumor}")
      crep <- "cancer_report/cancer_report_tables"
      regexes <- tibble::tribble(
        ~regex                                                                                           , ~fun                         ,
        glue("{crep}/{pref}-qc_summary\\.tsv\\.gz$")                                                     , "read_qcSum"                 ,
        glue("{crep}/purple/{pref}-purple_cnv_som_gene\\.tsv\\.gz$")                                     , "DOWNLOAD_ONLY-purplegene"   ,
        glue("{crep}/purple/{pref}-purple_cnv_som\\.tsv\\.gz$")                                          , "DOWNLOAD_ONLY-purplesom"    ,
        glue("{crep}/hrd/{pref}-hrdetect\\.tsv\\.gz$")                                                   , "read_hrdHrdetect"           ,
        glue("{crep}/sigs/{pref}-snv_2015\\.tsv\\.gz$")                                                  , "read_sigsTsv"               ,
        glue("{crep}/sigs/{pref}-snv_2020\\.tsv\\.gz$")                                                  , "read_sigsTsv"               ,
        glue("{crep}/sigs/{pref}-dbs\\.tsv\\.gz$")                                                       , "read_sigsTsv"               ,
        glue("{crep}/sigs/{pref}-indel\\.tsv\\.gz$")                                                     , "read_sigsTsv"               ,
        glue("smlv_somatic/filter/{libid_tumor}\\.pass\\.vcf\\.gz$")                                     , "DOWNLOAD_ONLY-smlvfiltvcf"  ,
        glue("smlv_somatic/filter/{libid_tumor}\\.pass\\.vcf\\.gz\\.tbi$")                               , "DOWNLOAD_ONLY-smlvfiltvcfi" ,
        glue("smlv_somatic/report/pcgr/{libid_tumor}\\.pcgr_acmg\\.grch38\\.json\\.gz$")                 , "read_pcgrJson"              ,
        glue("smlv_somatic/report/{libid_tumor}\\.somatic\\.variant_counts_process\\.json$")             , "read_smlvSomCounts"         ,
        # glue("smlv_somatic/report/pcgr/{libid_tumor}\\.pcgr_acmg\\.grch38\\.vcf\\.gz$")                  , "DOWNLOAD_ONLY-pcgrvcf"      ,
        # glue("smlv_somatic/report/pcgr/{libid_tumor}\\.pcgr_acmg\\.grch38\\.vcf\\.gz\\.tbi$")            , "DOWNLOAD_ONLY-pcgrvcfi"     ,
        glue("smlv_somatic/report/pcgr/{libid_tumor}\\.pcgr_acmg\\.grch38\\.snvs_indels\\.tiers\\.tsv$") , "DOWNLOAD_ONLY-pcgrtiers"    ,
        # glue("smlv_germline/report/cpsr/{libid_normal}\\.cpsr\\.grch38\\.vcf\\.gz$")                     , "DOWNLOAD_ONLY-cpsrvcf"      ,
        # glue("smlv_germline/report/cpsr/{libid_normal}\\.cpsr\\.grch38\\.vcf\\.gz\\.tbi$")               , "DOWNLOAD_ONLY-cpsrvcfi"     ,
        glue("sv_somatic/prioritise/{libid_tumor}\\.sv\\.prioritised\\.vcf\\.gz$")                       , "DOWNLOAD_ONLY-svpriovcf"    ,
        glue("sv_somatic/prioritise/{libid_tumor}\\.sv\\.prioritised\\.vcf\\.gz\\.tbi$")                 , "DOWNLOAD_ONLY-svpriovcfi"   ,
      )
      super$initialize(path = path, wname = wname, regexes = regexes)
      self$libid_tumor <- libid_tumor
      self$libid_normal <- libid_normal
    },
    #' @description Print details about the Workflow.
    #' @param ... (ignored).
    print = function(...) {
      res <- tibble::tribble(
        ~var           , ~value                               ,
        "path"         , private$.path                        ,
        "wname"        , private$.wname                       ,
        "filesystem"   , private$.filesystem                  ,
        "nregexes"     , as.character(nrow(private$.regexes)) ,
        "libid_tumor"  , self$libid_tumor                     ,
        "libid_normal" , self$libid_normal
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
      colnames(dat) <- tolower(colnames(dat))
      colnames(dat) <- gsub("\\.", "_", colnames(dat))
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
        dplyr::mutate(
          variable = tolower(.data$variable),
          variable = gsub("\\(|\\)", "", .data$variable),
          variable = gsub(" |-", "_", .data$variable)
        ) |>
        tidyr::pivot_wider(names_from = "variable", values_from = "value")
      tibble::tibble(name = "qcsum", data = list(dat[]))
    }
  ) # end public
)

#' sash Download Tidy and Write
#'
#' Downloads files from the `sash` workflow and writes them in a tidy format.
#'
#' @param path Path to directory with raw workflow results (from S3 or
#' local filesystem).
#' @param libid_tumor The LibraryID of the tumor sample.
#' @param libid_normal The LibraryID of the normal sample.
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
#' libid_tumor <- "PRJ230004"
#' p1_gds <- glue("gds://production/analysis_data/umccrise")
#' p <- file.path(p1_gds, "20240830ec648f40/L2300064__L2300063")
#' outdir <- file.path(sub("gds:/", "~/icav1/g", p))
#' d <- Wf_sash_download_tidy_write(
#'   path = p, libid_tumor = libid_tumor,
#'   outdir = outdir,
#'   dryrun = F
#' )
#' }
#' @export
Wf_sash_download_tidy_write <- function(
  path,
  libid_tumor,
  libid_normal,
  outdir,
  format = "rds",
  max_files = 1000,
  regexes = NULL,
  dryrun = FALSE
) {
  s <- Wf_sash$new(
    path = path,
    libid_tumor = libid_tumor,
    libid_normal = libid_normal
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
      prefix = glue("__{libid_tumor}"),
      format = format
    )
    return(d_write)
  }
  return(d_dl)
}
