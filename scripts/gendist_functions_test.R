# Dependencies:
library(tidyverse)
library(ggplot2)
library(purrr)
library(VennDiagram)

# Function to count aasubs per sequence ########################################
# Before this documentation is offically fomatted for as CRAN package, we can
# preview documentation using the docstring package:
# library(docstring)
# ?gendist_avg or docstring(gendist_avg)

gendist_avg <- function(test, focal, plot = FALSE) {
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
  #'   uses the "aaSubstitutions", "clade", and "seqName" (id) columns.
  #' @param plot an optional logical parameter to return a venn diagram gList object of
  #'   the average total aaSubstitutions. Object is created using draw.pairwise.venn() from
  #'   package "VennDiagram". Default is `FALSE`.
  #' @return returns a 13x13 dataframe object containing the breakdown of the
  #'   average comparison of test and focal, per genome segment. Or if
  #'   plot=TRUE, returns a list of 2: df = a 13x13 dataframe, plot = a gList
  #'   object.
  #' @export ??? this may be to identify the function when I put it in the
  #'   correct format? But the keyboard shortcut for template geneation doesnt
  #'   put it at the bottom..
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
  #' "clade" = c("21K"),
  #' "aaSubstitutions" = c("E:XXX,E:XX1,E:XX3")
  #' )
  #'
  #' # Execute function---------------------------------------------------------
  #' results <- gendist_avg(test, focal)
  #' results
  #'
  #' results_with_plot <- gendist_avg(test, focal, plot=TRUE)
  #' view(results_with_plot$df)
  #' grid.draw(results_with_plot$plot)

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
      dplyr::mutate(!!i := str_count(aaSubstitutions, paste0(i, ":")))
  }

  # total number of aasubs
  subs <- subs %>%
    dplyr::mutate(total = rowSums(dplyr::select(., E:S), na.rm = TRUE))

  # get mean column
  mean_subs <- subs %>%
    dplyr::select(-c(aaSubstitutions)) %>%
    summarise_all(mean) %>%
    t() %>%
    as.data.frame() %>%
    rename(testdf_mean = V1)

  # get sd column
  sd_subs <- subs %>%
    dplyr::select(-c(aaSubstitutions)) %>%
    summarise_all(sd) %>%
    t() %>%
    as.data.frame() %>%
    rename(testdf_sd = V1)

  # combine
  mean_subs <- dplyr::bind_cols(mean_subs, sd_subs)

  # count focal set ------------------------------------------------------------

  # trim data
  focal <- focal %>%
    dplyr::select(seqName, clade, aaSubstitutions)

  # counts number of aasubs per segment
  for (i in segments) {
    focal <- focal %>%
      dplyr::mutate(!!i := str_count(aaSubstitutions, paste0(i, ":")))
  }

  # total number of aasubs
  focal <- focal %>%
    dplyr::mutate(total = rowSums(dplyr::select(., E:S), na.rm = TRUE))

  # get mean column
  focal_mean <- focal %>%
    dplyr::select(-c(aaSubstitutions, clade, seqName)) %>%
    summarise_all(mean) %>%
    t() %>%
    as.data.frame() %>%
    rename(focal_mean = V1)

  # get sd column
  focal_sd <- focal %>%
    dplyr::select(-c(aaSubstitutions, clade, seqName)) %>%
    summarise_all(sd) %>%
    t() %>%
    as.data.frame() %>%
    rename(focal_sd = V1)

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
    dplyr::mutate(across(.cols = c(16:last_col()), ~ purrr::map2(test_exp, ., ~ gsub(.x, "", .y))))

  strains <- focal$seqName # this will iterate for every unique seqName

  df <- data.frame() # empty df to store results

  for (i in strains) { # iterate over every row, uniquely identified by strain
    iter <- subs_string %>%
      dplyr::select(all_of(i))

    for (j in segments) { # iterate over every segment for every strain
      iter <- iter %>%
        dplyr::mutate(!!j := str_count(!!sym(i), paste0(j, ":")))
    }

    # Sum total unique aasubs
    iter <- iter %>%
      dplyr::mutate(total = rowSums(dplyr::select(., E:S), na.rm = TRUE))

    iter <- iter %>%
      dplyr::select(all_of(segments), total)

    df <- bind_rows(df, iter) # bind results
  }

  # get mean column
  focal_rawdist_mean <- df %>% # get mean of all pairwise comparisons
    summarize_all(mean) %>%
    t() %>%
    as.data.frame() %>%
    rename(focal_rawdist_mean = V1)

  # get sd column
  focal_rawdist_sd <- df %>%
    summarize_all(sd) %>%
    t() %>%
    as.data.frame() %>%
    rename(focal_rawdist_sd = V1)

  # combine
  focal_unique <- dplyr::bind_cols(focal_rawdist_mean, focal_rawdist_sd)

  # unique to test set ------------------------------------------------------------

  # remove focal patterns from test strings
  subs_string <- dplyr::bind_cols(subs, focal_patterns) %>%
    dplyr::mutate(across(.cols = c(16:last_col()), ~ purrr::map2(., aaSubstitutions, ~ gsub(.x, "", .y))))

  strains <- focal$seqName

  df <- data.frame() # empty df to store results

  for (i in strains) { # iterate over every row, uniquely identified by seqName
    iter <- subs_string %>%
      dplyr::select(all_of(i))

    for (j in segments) { # iterate over every segment for every seqName
      iter <- iter %>%
        dplyr::mutate(!!j := str_count(!!sym(i), paste0(j, ":")))
    }

    # Sum total unique aasubs
    iter <- iter %>%
      dplyr::mutate(total = rowSums(dplyr::select(., E:S), na.rm = TRUE))

    iter <- iter %>%
      dplyr::select(all_of(segments), total)

    df <- bind_rows(df, iter) # bind results
  }

  # get mean column
  testdf_rawdist_mean <- df %>% # get mean of all pairwise comparisons
    summarize_all(mean) %>%
    t() %>%
    as.data.frame() %>%
    rename(testdf_rawdist_mean = V1)

  # get sd column
  testdf_rawdist_sd <- df %>%
    summarize_all(sd) %>%
    t() %>%
    as.data.frame() %>%
    rename(testdf_rawdist_sd = V1)

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

  avg_dist_per_seq <- rownames_to_column(avg_dist_per_seq, "segment")

  if (plot == FALSE) {
    return(avg_dist_per_seq)
  }

  if (plot == TRUE) {
    vennplot <- draw.pairwise.venn(
      area1 = avg_dist_per_seq$testdf_mean[avg_dist_per_seq$segment == "total"], # Create pairwise venn diagram
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

# Function to parse aasubs per sequence ########################################
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

  if(is.null(focal)){

    subs<- subs %>%
      dplyr::select(-aaSubstitutions)

    return(subs)

  }

  if(!is.null(focal)){

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

# Function to count overall aasubs in sample (not per sequence) ###############

gendist_all <- function(data, focal, nc) {
  # create vector of genome segments
  segments <- c("E", "M", "N", "ORF1a", "ORF1b", "ORF3a", "ORF6", "ORF7a", "ORF7b", "ORF8", "ORF9b", "S")

  # trim down data
  overall <- data %>%
    dplyr::select(aaSubstitutions) %>%
    dplyr::filter(!is.na(aaSubstitutions))

  # test set -------------------------------------------------------------------

  # pull out unique mutations in test sample
  test <- unlist(str_split(overall$aaSubstitutions, ","))
  unique_test <- unique(test)
  unique_test <- paste0(unique_test, collapse = ",")
  overall <- data.frame(unique_test)
  overall_test <- overall

  # count number of aasubs per genome segment
  for (i in segments) {
    overall_test <- overall_test %>%
      dplyr::mutate(!!i := str_count(unique_test, paste0(i, ":")))
  }

  # Sum counts and create column
  overall_test <- overall_test %>%
    dplyr::mutate(total = rowSums(dplyr::select(., E:S), na.rm = TRUE)) %>%
    dplyr::select(-unique_test) %>%
    t() %>%
    as.data.frame() %>%
    rename(testdf_rel_wt = V1)

  # focal set --------------------------------------------------------------

  # trim data
  focal <- focal %>%
    dplyr::select(seqName, clade, aaSubstitutions)

  # count number of aasubs per genome segment
  for (i in segments) {
    focal <- focal %>%
      dplyr::mutate(!!i := str_count(aaSubstitutions, paste0(i, ":")))
  }

  # total number
  focal <- focal %>%
    dplyr::mutate(total = rowSums(dplyr::select(., E:S), na.rm = TRUE))

  # get focal count column
  focal_counts <- focal %>%
    dplyr::filter(clade == nc) %>% # for BA.1
    dplyr::select(-c(seqName, clade, aaSubstitutions)) %>%
    t() %>%
    as.data.frame() %>%
    rename(focal_rel_wt = V1)

  # create focal comparison string
  lin <- focal %>%
    dplyr::filter(clade == nc) %>%
    dplyr::select(aaSubstitutions)

  focal_string <- lin[[1]]
  focal_vec <- unlist(str_split(focal_string, ","))
  focal_pattern <- paste0(focal_vec, collapse = "|")

  # unique test set ----------------------------------------------------------

  # trim down data
  subs_unique <- overall %>%
    dplyr::mutate(aasubs_unique = gsub(focal_pattern, "", unique_test))

  # count number of unique aasubs per segment
  for (i in segments) {
    subs_unique <- subs_unique %>%
      dplyr::mutate(!!i := str_count(aasubs_unique, paste0(i, ":")))
  }

  # sum aasubs and create column
  overall_unique <- subs_unique %>%
    dplyr::mutate(total = rowSums(dplyr::select(., E:S), na.rm = TRUE)) %>%
    dplyr::select(-c(unique_test, aasubs_unique)) %>%
    t() %>%
    as.data.frame() %>%
    rename(testdf_unique = V1)

  # create table -----------------------------------------------------------------
  overall_dist <- dplyr::bind_cols(overall_test, overall_unique, focal_counts) %>%
    dplyr::mutate(testdf_focal_shared = testdf_rel_wt - testdf_unique) %>%
    dplyr::mutate(focal_unique = focal_rel_wt - testdf_focal_shared) %>%
    dplyr::select(testdf_rel_wt, testdf_unique, testdf_focal_shared, focal_unique, focal_rel_wt)

  overall_dist <- round(overall_dist, 1)

  overall_dist <- rownames_to_column(overall_dist, "segment")

  return(overall_dist)
}


# Function to parse all aasubs found in sample (overall, not per sequence)######

gendist_parse_all <- function(data, focal, nc) {
  # create vector of genome segments
  segments <- c("E", "M", "N", "ORF1a", "ORF1b", "ORF3a", "ORF6", "ORF7a", "ORF7b", "ORF8", "ORF9b", "S")

  # trim down data
  overall <- data %>%
    dplyr::select(aaSubstitutions) %>%
    dplyr::filter(!is.na(aaSubstitutions))

  # test set -------------------------------------------------------------------

  # pull out unique aasubs in test sample
  test <- unlist(str_split(overall$aaSubstitutions, ","))
  unique_test <- unique(test)
  unique_test <- paste0(unique_test, collapse = ",")
  overall <- data.frame(unique_test)
  overall_test <- overall

  # create focal comparison string
  lin <- focal %>%
    dplyr::filter(clade == nc) %>%
    dplyr::select(aaSubstitutions)

  focal_string <- lin[[1]]
  focal_vec <- unlist(str_split(focal_string, ","))

  # parse aasubs per segment
  for (i in segments) {
    overall_test <- overall_test %>%
      dplyr::mutate(!!i := stringr::str_extract_all(unique_test, paste0(i, ":.*?(?=,)")))
  }

  # parse unique aasubs per segment (not shared with focal)
  overall_test$unique_E <- lapply(overall_test$E, function(x) setdiff(unlist(x), focal_vec))
  overall_test$unique_M <- lapply(overall_test$M, function(x) setdiff(unlist(x), focal_vec))
  overall_test$unique_N <- lapply(overall_test$N, function(x) setdiff(unlist(x), focal_vec))
  overall_test$unique_ORF1a <- lapply(overall_test$ORF1a, function(x) setdiff(unlist(x), focal_vec))
  overall_test$unique_ORF1b <- lapply(overall_test$ORF1b, function(x) setdiff(unlist(x), focal_vec))
  overall_test$unique_ORF3a <- lapply(overall_test$ORF3a, function(x) setdiff(unlist(x), focal_vec))
  overall_test$unique_ORF6 <- lapply(overall_test$ORF6, function(x) setdiff(unlist(x), focal_vec))
  overall_test$unique_ORF7a <- lapply(overall_test$ORF7a, function(x) setdiff(unlist(x), focal_vec))
  overall_test$unique_ORF7b <- lapply(overall_test$ORF7b, function(x) setdiff(unlist(x), focal_vec))
  overall_test$unique_ORF8 <- lapply(overall_test$ORF8, function(x) setdiff(unlist(x), focal_vec))
  overall_test$unique_ORF9b <- lapply(overall_test$ORF9b, function(x) setdiff(unlist(x), focal_vec))
  overall_test$unique_S <- lapply(overall_test$S, function(x) setdiff(unlist(x), focal_vec))
  return(overall_test)
}


# Venn Diagrams ----------------------------------------------------------------
venn_dist <- function(data, focal, nc) {
  # create vector of genome segments
  segments <- c("E", "M", "N", "ORF1a", "ORF1b", "ORF3a", "ORF6", "ORF7a", "ORF7b", "ORF8", "ORF9b", "S")

  # trim down data
  overall <- data %>%
    dplyr::select(aaSubstitutions) %>%
    dplyr::filter(!is.na(aaSubstitutions))

  # pull out unique mutations in test sample
  test <- unlist(str_split(overall$aaSubstitutions, ","))
  unique_test <- unique(test)

  # trim data
  focal <- focal %>%
    dplyr::select(seqName, clade, aaSubstitutions)

  # create focal comparison string
  lin <- focal %>%
    dplyr::filter(clade == nc) %>%
    dplyr::select(aaSubstitutions)

  focal_string <- lin[[1]]
  focal_vec <- unlist(str_split(focal_string, ","))

  venn_diagrams <- list()
  # per genome segment
  for (i in segments) {
    # parse
    focal2 <- focal_vec[grepl(paste0(i, ":"), focal_vec)]
    test2 <- unique_test[grepl(paste0(i, ":"), unique_test)]
    ventest <- list(test2, focal2)
    # plot
    p <- ggVennDiagram(ventest,
      category.names = c("Test", "Focal")
    ) +
      ggtitle(paste0("Genome Segment: ", i)) +
      guides(fill = guide_legend(title = "aasubs"))
    venn_diagrams[[i]] <- p
  }
  # overall
  ventest <- list(unique_test, focal_vec)
  p <- ggVennDiagram(ventest,
    category.names = c("Test", "Focal")
  ) +
    ggtitle("Overall Genome") +
    guides(fill = guide_legend(title = "aasubs"))
  venn_diagrams <- c(venn_diagrams, list(overall = p))
  return(venn_diagrams)
}
