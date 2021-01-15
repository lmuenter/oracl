###########################
#
# Utility Functions
#
###########################

# purpose: Utility Functions for ORA testing
# date: 15.01.2021
# author: Lukas Muenter
# note: only tested for Arabidopsis thaliana!

#' Add list names as column and combine
#'
#' Takes a list of ORAs, adds column "grouping" with element names and combines to a single ORA dataframe
#' @param als A list of ORA dataframes
#' @return One ORA dataframe
#' @export

oracl_list_to_df <- function(als){

  ls.names = names(als)
  df = lapply(

    ls.names,
    function(x,y){
      x[[y]] %>%
        mutate("grouping" = y)
    },
    x = als

  ) %>%
    do.call("rbind", .)

}

