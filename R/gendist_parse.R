#' Parse Amino Acid Substitutions Per Genome Segment
#'
#' #' `gendist_parse()` calculates parse the aaSubs using the "aaSubstitutions"
#' column generated from a Nextclade metadata output
#' (<https://clades.nextstrain.org/>). Output includes overall, unique_test, and
#' unique_focal.
#'
#' @param test a dataframe containing nextclade metadata output for all
#'   sequences in the test set (often circulating lineages within a certain
#'   location and timeframe). Function works on the aaSubstitutions column.
#'
#' @param focal optional parameter that takes a dataframe containing nextclade
#'   metadata output for a single focal sequence.
#'
#' @return if no focal is supplied, the function returns a dataframe of amino
#'   acid substitutions parsed by genome segment for each sequence in test (13
#'   columns total). If a focal parameter is supplied, the function returns an
#'   expanded dataframe with two additional parses: unique to test, and unique
#'   to focal (37 columns total).
#'
#' @export
#'
#' @examples
#'
#' # Create test df-----------------------------------------------------------
#' test <- data.frame(
#' "aaSubstitutions" = c("E:XXX,E:XX1,S:XXX", "E:XXX,E:XX2,M:XXX"),
#' "seqName" = c(1, 2)
#' )
#'
#' # Create focal df----------------------------------------------------------
#' focal <- data.frame(
#' "aaSubstitutions" = c("E:XXX,E:XX1,E:XX3,S:XXX,MXX1"),
#' "seqName" = c("3")
#' )
#'
#' # Execute function---------------------------------------------------------
#' results <- gendist_parse(test, focal)
#' results
#'
#' results_with_focal <- gendist_parse(test, focal)
#' results_with_focal
#'
gendist_parse <- function(test, focal=NULL) {
  # create vector of genome segments
  segments <- c("E", "M", "N", "ORF1a", "ORF1b", "ORF3a", "ORF6", "ORF7a", "ORF7b", "ORF8", "ORF9b", "S")

  # trim down test
  subs <- test %>%
    dplyr::select(seqName, aaSubstitutions) %>%
    dplyr::filter(!is.na(aaSubstitutions))

  for (i in segments) {
    subs<- subs %>%
      dplyr::mutate(!!i := stringr::str_extract_all(aaSubstitutions, paste0(i, ":.*?(?=,|$)")))
  }

  if(is.null(focal)){ # return only the test parse if no focal is supplied

    subs<- subs %>%
      dplyr::select(-aaSubstitutions)

    return(subs)

  }

  if(!is.null(focal)){ # include a parse for unique to test and unique to focal

    # get string and regex patterns
    focal <- focal %>%
      dplyr::mutate(
        focal_string = aaSubstitutions,
        focal_pattern = gsub(",", "|", aaSubstitutions)
      )
    # focal pattern
    focal_patterns <- focal %>%
      dplyr::select(seqName, focal_pattern) %>%
      dplyr::mutate(seqName = "focal_exp") %>%
      tidyr::pivot_wider(
        names_from = seqName,
        values_from = focal_pattern
      )

    # focal string
    focal_strings <- focal %>%
      dplyr::select(seqName, focal_string) %>%
      dplyr::mutate(seqName = "focal") %>%
      tidyr::pivot_wider(
        names_from = seqName,
        values_from = focal_string
      )

    # get test regex pattern
    subs <- subs %>%
      dplyr::mutate(test_exp = gsub(",", "|", aaSubstitutions))

    # remove test patterns from focal strings
    subs_string <- dplyr::bind_cols(subs, focal_strings)%>%
      dplyr::mutate(focal_unique = purrr::map2(test_exp, focal, ~ gsub(.x, "", .y)))

    # remove focal pattern from test strings
    subs_string <- dplyr::bind_cols(subs_string, focal_patterns)%>%
      dplyr::mutate(test_unique = purrr::map2(focal_exp, aaSubstitutions, ~ gsub(.x, "", .y)))

    # Parse unique aasubs per genome segment in test
    for (i in segments) {
      subs_string <- subs_string %>%
        dplyr::mutate(!!paste0(i, "_test_unique") := stringr::str_extract_all(test_unique, paste0(i, ":.*?(?=,|$)")))
    }

    # Parse unique aasubs per genome segment in focal
    for (i in segments) {
      subs_string <- subs_string %>%
        dplyr::mutate(!!paste0(i, "_focal_unique") := stringr::str_extract_all(focal_unique, paste0(i, ":.*?(?=,|$)")))
    }

    subs_string <- subs_string %>%
      dplyr::select(-c(aaSubstitutions, test_exp, focal, focal_unique, focal_exp, test_unique))

    return(subs_string)
  }

}
