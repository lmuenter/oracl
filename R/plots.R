############################
#
# Fancy Plots
#
############################

# purpose: adds fancy plots to visualise enriched GO-Terms
# date: 27.07.2021
# author: Lukas Muenter

#' Volcano Plot
#'
#' Creates a volcano plot from the output of `oracl::oraclient`
#' @param df A dataframe generated using `oracl::oraclient`
#' @param x NUMERIC column to be used for the x axis of the plot. Will be log2-transformed. (default: `"fold_change"`)
#' @param y Significance metric to be used for the x axis of the plot. Will be -log-transformed. Options: `c("p.value", "fdr")` (default: `"fdr"`)
#' @param top_n Number of terms that should be labelled. Default: 5.
#' @param trunc To which size should Term-label be truncated? Defaults to `20`
#' @param side From which side should truncation happen? Defauls to `left`
#' @param include_ID should GO_ID be included in plot labels? Defaults to `TRUE`
#' @import ggplot2
#' @importFrom tibble add_column
#' @importFrom ggrepel geom_text_repel
#' @export
volcanoracl = function(df, x = "fold_change", y = "fdr", top_n = 5, trunc = 20, side = "left", include_ID = TRUE){

  ## check dataframe
  if(!is.data.frame(df)) {

    stop("Error: List was provided as input. input should be a dataframe. Consider using oracl::oracl_list_to_df prior to calling this function." )

  }

  ## check input for Y
  if(length(y) > 1 | length(y) < 1){

    stop("Error: More than one element provided. `y` should be either 'fdr' or 'p.value'")

  } else if(is.character(y) == FALSE){

    stop("Error: Input is not a character. `y` should be either 'fdr' or 'p.value'")

  } else if(!y %in% c("p.value", "fdr")){

    stop("Error: `y` should be either 'fdr' or 'p.value'")

  }

  ## check input for X
  if(length(x) > 1 | length(x) < 1){

    stop("Error: More than one element provided. `x` should specify ONE numeric column.")

  } else if(is.numeric(df[[x]]) == FALSE){

    stop("Error: `x` should specify one numeric column.")

  }

  ## truncate labels, if desired, include GO_ID
  if(include_ID){

    df <- df %>%
      mutate(
        "label" = paste(str_trunc(label, 20, side = side),
                        " (",
                        GO_ID,
                        ")",
                        sep = "")
      )
  } else {

    df <- df %>%
      mutate(
        "label" = str_trunc(label, 20, side = side)
      )

  }


  ## if column "grouping" doesn't exist, add the group
  cols <- c(grouping = NA_real_)
  df <- add_column(df, !!!cols[setdiff(names(cols), names(df))])

  ## Get top terms per group
  top_terms <- df %>%
    group_by(grouping) %>%
    slice_max(-log(fdr), n = top_n) %>%
    mutate("labelthis" = label) %>%
    select(label, labelthis)

  ## Plot
  p <- df %>%
    left_join(top_terms) %>%
    mutate(labelthis = ifelse(!is.na(labelthis), labelthis, "")) %>%
    ggplot(aes(x = log2(fold_change), y = -log(fdr), label = labelthis, color = grouping)) +
    geom_point(show.legend = FALSE) +
    geom_text_repel(color = "black") +
    theme_bw()

  ## plot the darn thing
  p

}

#' Dot Plot
#'
#' Creates a dot plot from the output of `oracl::oraclient`.
#' @param df A dataframe generated using `oracl::oraclient`.
#' @param top_n if set, number of terms that should be labelled. Default: `50`. Select `NULL` for no filter at all (not recommended, will get cluttered).
#' @import ggplot2
#' @param trunc To which size should Term-label be truncated? Defaults to `20`
#' @param side From which side should truncation happen? Defauls to `left`
#' @param include_ID should GO_ID be included in plot labels? Defaults to `TRUE`
#' @importFrom tibble add_column
#' @importFrom stringr str_trunc
#' @export
oraclot = function(df, top_n = 50, trunc = 20, side = "left", include_ID = TRUE){

  ## check dataframe
  if(!is.data.frame(df)) {

    stop("Error: List was provided as input. input should be a dataframe. Consider using oracl::oracl_list_to_df prior to calling this function." )

  }

  ## truncate labels, if desired, include GO_ID
  if(include_ID){

    df <- df %>%
      mutate(
        "label" = paste(str_trunc(label, 20, side = side),
                        " (",
                        GO_ID,
                        ")",
                        sep = "")
      )
  } else {

    df <- df %>%
      mutate(
        "label" = str_trunc(label, 20, side = side)
      )

  }

  ## get label order
  label_order <- df %>% arrange(fold_change) %>% pull(label)

  ## if column "grouping" doesn't exist, add the group
  cols <- c(grouping = NA_real_)
  df <- add_column(df, !!!cols[setdiff(names(cols), names(df))])

  ## Get top terms per group
  top_terms <- df %>%
    group_by(grouping) %>%
    slice_max(-log(fdr), n = top_n) %>%
    mutate("labelthis" = label) %>%
    select(label, labelthis)

  ## create basic plot.
  p = df %>%
    left_join(top_terms) %>%
    mutate(labelthis = ifelse(!is.na(labelthis), labelthis, "")) %>%
    filter(labelthis != "") %>%
    mutate("labelthis" = factor(labelthis, levels = label_order)) %>%
    ggplot(aes(x = log2(fold_change), y = labelthis, size = -log(p.value), color = -log(fdr))) +
    geom_point() +
    theme_bw()

  ## plot the darn thing
  p

}


