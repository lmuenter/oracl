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
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @export
volcanoracl = function(df, x = "fold_change", y = "fdr", top_n = 5){

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

  ## get topN genes
  top_genes <- df %>%
    slice_max(-log(fdr), n = top_n) %>%
    pull(label)

  ## create basic plot
  p <- df %>%
    mutate(labelthis = ifelse(label %in% top_genes, label, ""),
           grouping = ifelse(!"grouping" %in% colnames(.), "a", grouping)) %>%
    ggplot(aes(x = log2(fold_change), y = -log(fdr), label = labelthis, color = grouping)) +
    geom_point(show.legend = FALSE) +
    geom_text_repel(color = "black") +
    theme_bw()

  ## if grouped dataframe, facet the plot
  if("grouping" %in% colnames(df)){

    p <- p + facet_wrap(grouping ~ .)

  }

  ## plot the darn thing
  p

}

#' Dot Plot
#'
#' Creates a dot plot from the output of `oracl::oraclient`.
#' @param df A dataframe generated using `oracl::oraclient`.
#' @param top_n if set, number of terms that should be labelled. Default: `50`. Select `NULL` for no filter at all (not recommended, will get cluttered).
#' @import ggplot2
#' @export
oraclot = function(df, top_n = 50){

  ## check dataframe
  if(!is.data.frame(df)) {

    stop("Error: List was provided as input. input should be a dataframe. Consider using oracl::oracl_list_to_df prior to calling this function." )

  }

  ## get topN genes
  if(!is.null(top_n)){

    top_terms <- df %>%
      slice_max(fold_change, n = top_n) %>%
      pull(label)

  } else {

    top_terms <- NULL

  }

  ## get label order
  label_order <- df %>% arrange(fold_change) %>% pull(label)

  ## create basic plot. If top_terms are provided, only use those
  if(is.null(top_terms)){

    p <- df %>%
      mutate("label" = factor(label, levels = label_order)) %>%
      ggplot(aes(x = log2(fold_change), y = label, size = -log(p.value), color = -log(fdr))) +
      geom_point() +
      theme_bw()

  } else {

    p <- df %>%
      mutate("label" = factor(label, levels = label_order)) %>%
      filter(label %in% top_terms) %>%
      ggplot(aes(x = log2(fold_change), y = label, size = -log(p.value), color = -log(fdr))) +
      geom_point() +
      theme_bw()

  }

  ## if grouped dataframe, facet the plot
  if("grouping" %in% colnames(df)){

    p <- p + facet_wrap(grouping ~ .)

  }

  ## plot the darn thing
  p

}


