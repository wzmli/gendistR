#' Calculate mean genetic distance
#'
#' `gendist_avg()` calculates the average of various genetic distance metrics
#' for two groups of SARS-CoV-2 sequences using the "aaSubstitutions" column
#' generated from a Nextclade metadata output
#' (<https://clades.nextstrain.org/>).
#'
#' For a detailed explanation of each metric produced, see the vignette.
#'
#' @param test a dataframe containing nextclade metadata output for all
#'   sequences in the test set (often circulating lineages within a certain
#'   location and timeframe). Function works on the aaSubstitutions column.
#' @param focal a dataframe containing nextclade metadata output for all focal
#'   sequences (can be one emerging lineage, or a set of multiple). Function
#'   uses the "aaSubstitutions", and "seqName" (id) columns.
#' @param plot an optional logical parameter to return a venn diagram gList object of
#'   the average total aaSubstitutions. Object is created using VennDiagram::draw.pairwise.venn() from
#'   package "VennDiagram". Default is `FALSE`.
#' @return returns a 13x13 dataframe object containing the breakdown of the
#'   average comparison of test and focal, per genome segment. Or if
#'   plot=TRUE, returns a list of 2: df = a 13x13 dataframe, plot = a gList
#'   object.
#' @importFrom magrittr %>%
#' @importFrom rlang :=
#' @importFrom rlang !!
#' @importFrom stats sd
#' @export
#' @examples
#' # Simple example focusing only on the E segment----------------------------
#'
#' # Create test df-----------------------------------------------------------
#' test <- data.frame(
#' "aaSubstitutions" = c("E:XXX,E:XX1", "E:XXX,E:XX2")
#' )
#'
#' # Create focal df----------------------------------------------------------
#' focal <- data.frame(
#' "seqName" = c("1"),
#' "aaSubstitutions" = c("E:XXX,E:XX1,E:XX3")
#' )
#'
#' # Execute function---------------------------------------------------------
#' results <- gendist_avg(test, focal)
#' results
#'
#' results_with_plot <- gendist_avg(test, focal, plot=TRUE)
#' tibble::view(results_with_plot$df)
#' grid::grid.draw(results_with_plot$plot)
#'
gendist_avg <- function(test, focal, plot = FALSE) {
  # create vector of genome segments
  segments <- c("E", "M", "N", "ORF1a", "ORF1b", "ORF3a", "ORF6", "ORF7a", "ORF7b", "ORF8", "ORF9b", "S")

  # trim down data
  subs <- test %>%
    dplyr::select(aaSubstitutions) %>%
    dplyr::filter(!is.na(aaSubstitutions))

  # count test set -------------------------------------------------------------

  # count number of aasubs per genome segment
  for (i in segments) {
    subs <- subs %>%
      dplyr::mutate(!!i := stringr::str_count(aaSubstitutions, paste0(i, ":")))
  }

  # total number of aasubs
  subs <- subs %>%
    dplyr::mutate(total = rowSums(dplyr::select(., E:S), na.rm = TRUE))

  # get mean column
  mean_subs <- subs %>%
    dplyr::select(-c(aaSubstitutions)) %>%
    dplyr::summarize_all(mean) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::rename(testdf_mean = V1)

  # get sd column
  sd_subs <- subs %>%
    dplyr::select(-c(aaSubstitutions)) %>%
    dplyr::summarize_all(sd) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::rename(testdf_sd = V1)

  # combine
  mean_subs <- dplyr::bind_cols(mean_subs, sd_subs)

  # count focal set ------------------------------------------------------------

  # trim data
  focal <- focal %>%
    dplyr::select(seqName, aaSubstitutions)

  # counts number of aasubs per segment
  for (i in segments) {
    focal <- focal %>%
      dplyr::mutate(!!i := stringr::str_count(aaSubstitutions, paste0(i, ":")))
  }

  # total number of aasubs
  focal <- focal %>%
    dplyr::mutate(total = rowSums(dplyr::select(., E:S), na.rm = TRUE))

  # get mean column
  focal_mean <- focal %>%
    dplyr::select(-c(aaSubstitutions, seqName)) %>%
    dplyr::summarize_all(mean) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::rename(focal_mean = V1)

  # get sd column
  focal_sd <- focal %>%
    dplyr::select(-c(aaSubstitutions, seqName)) %>%
    dplyr::summarize_all(sd) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::rename(focal_sd = V1)

  # combine
  mean_focal <- dplyr::bind_cols(focal_mean, focal_sd)

  # count unique focal ---------------------------------------------------------

  # get string and regex patterns
  focal <- focal %>%
    dplyr::mutate(
      focal_string = aaSubstitutions,
      focal_pattern = gsub(",", "|", aaSubstitutions)
    )
  # focal pattern
  focal_patterns <- focal %>%
    dplyr::select(seqName, focal_pattern) %>%
    tidyr::pivot_wider(
      names_from = seqName,
      values_from = focal_pattern
    )

  # focal string
  focal_strings <- focal %>%
    dplyr::select(seqName, focal_string) %>%
    tidyr::pivot_wider(
      names_from = seqName,
      values_from = focal_string
    )

  # get test regex pattern
  subs <- subs %>%
    dplyr::mutate(test_exp = gsub(",", "|", aaSubstitutions))

  # remove test patterns from focal strings
  subs_string <- dplyr::bind_cols(subs, focal_strings) %>%
    dplyr::mutate(dplyr::across(.cols = c(16:tidyr::last_col()), ~ purrr::map2(test_exp, ., ~ gsub(.x, "", .y))))

  strains <- focal$seqName # this will iterate for every unique seqName

  df <- data.frame() # empty df to store results

  for (i in strains) { # iterate over every row, uniquely identified by strain
    iter <- subs_string %>%
      dplyr::select(tidyr::all_of(i))

    for (j in segments) { # iterate over every segment for every strain
      iter <- iter %>%
        dplyr::mutate(!!j := stringr::str_count(!!rlang::sym(i), paste0(j, ":")))
    }

    # Sum total unique aasubs
    iter <- iter %>%
      dplyr::mutate(total = rowSums(dplyr::select(., E:S), na.rm = TRUE))

    iter <- iter %>%
      dplyr::select(tidyr::all_of(segments), total)

    df <- dplyr::bind_rows(df, iter) # bind results
  }

  # get mean column
  focal_rawdist_mean <- df %>% # get mean of all pairwise comparisons
    dplyr::summarize_all(mean) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::rename(focal_rawdist_mean = V1)

  # get sd column
  focal_rawdist_sd <- df %>%
    dplyr::summarize_all(sd) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::rename(focal_rawdist_sd = V1)

  # combine
  focal_unique <- dplyr::bind_cols(focal_rawdist_mean, focal_rawdist_sd)

  # unique to test set ------------------------------------------------------------

  # remove focal patterns from test strings
  subs_string <- dplyr::bind_cols(subs, focal_patterns) %>%
    dplyr::mutate(dplyr::across(.cols = c(16:tidyr::last_col()), ~ purrr::map2(., aaSubstitutions, ~ gsub(.x, "", .y))))

  strains <- focal$seqName

  df <- data.frame() # empty df to store results

  for (i in strains) { # iterate over every row, uniquely identified by seqName
    iter <- subs_string %>%
      dplyr::select(tidyr::all_of(i))

    for (j in segments) { # iterate over every segment for every seqName
      iter <- iter %>%
        dplyr::mutate(!!j := stringr::str_count(!!rlang::sym(i), paste0(j, ":")))
    }

    # Sum total unique aasubs
    iter <- iter %>%
      dplyr::mutate(total = rowSums(dplyr::select(., E:S), na.rm = TRUE))

    iter <- iter %>%
      dplyr::select(tidyr::all_of(segments), total)

    df <- dplyr::bind_rows(df, iter) # bind results
  }

  # get mean column
  testdf_rawdist_mean <- df %>% # get mean of all pairwise comparisons
    dplyr::summarize_all(mean) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::rename(testdf_rawdist_mean = V1)

  # get sd column
  testdf_rawdist_sd <- df %>%
    dplyr::summarize_all(sd) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::rename(testdf_rawdist_sd = V1)

  # combine
  test_unique <- dplyr::bind_cols(testdf_rawdist_mean, testdf_rawdist_sd)

  # create table -----------------------------------------------------------------

  avg_dist_per_seq <- dplyr::bind_cols(mean_subs, test_unique, focal_unique, mean_focal) %>%
    dplyr::mutate(testdf_focal_shared = testdf_mean - testdf_rawdist_mean) %>%
    dplyr::select(testdf_mean, testdf_sd, testdf_rawdist_mean, testdf_rawdist_sd, testdf_focal_shared, focal_rawdist_mean, focal_rawdist_sd, focal_mean, focal_sd) %>%
    dplyr::mutate(
      jaccard_mean = 1 - (testdf_focal_shared / (testdf_rawdist_mean + testdf_focal_shared + focal_rawdist_mean)),
      focal_propdist_mean = focal_rawdist_mean / focal_mean,
      testdf_propdist_mean = testdf_rawdist_mean / testdf_mean
    )

  avg_dist_per_seq <- round(avg_dist_per_seq, 2)

  avg_dist_per_seq <- tibble::rownames_to_column(avg_dist_per_seq, "segment")

  if (plot == FALSE) {
    return(avg_dist_per_seq)
  }

  if (plot == TRUE) {
    # Create pairwise venn diagram
    vennplot <- VennDiagram::draw.pairwise.venn(
      area1 = avg_dist_per_seq$testdf_mean[avg_dist_per_seq$segment == "total"],
      area2 = avg_dist_per_seq$focal_mean[avg_dist_per_seq$segment == "total"],
      cross.area = avg_dist_per_seq$testdf_focal_shared[avg_dist_per_seq$segment == "total"],
      category = c("Testdf", "Focal"),
      fill = c("red", "blue"),
      cat.pos = c(190, 170),
      print.mode = c("raw", "percent")
    )
    result <- list(df = avg_dist_per_seq, plot = vennplot)
  }
}
