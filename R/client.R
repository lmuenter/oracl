###########################
#
# PANTHER client
#
###########################

# purpose: simple client for sending a geneset and a background geneset to panther.org for an overrepresentation test
# date: 15.01.2021
# author: Lukas Muenter
# note: only tested for Arabidopsis thaliana!

#' Panther GO term enrichment
#'
#' A client for the panther webservice. Takes a list of genes and a background geneset and sends them to panther for a classic Overrepresentation Analysis (ORA)
#' @param x A list of gene identifiers (character vector)
#' @param bg A list of background genes (character vector)
#' @param taxon Name of organism to be studied (Defaults to: Athaliana (ID: 3702))
#' @param ontology Panther annotation data type of interest (Defaults to "GO:0008150" for BP, also possible: "GO:0003674" for MF, "GO:0005575" for CC)
#' @param enrichmentTestType Type of enrichment test (Defaults to: "FISHER", also possible: BINOMIAL)
#' @param correction FP-Correction method (Defaults to: "BONFERRONI", also possible: FDR)
#' @param fdr.thresh FDR cutoff (default: 0.1)
#' @param p.thresh P-value cutoff to be applied (default: p < 0.05)
#' @return A dataframe
#' @export

oracl <- function(x, bg, ontology = "bp", taxon = "Athaliana",  enrichmentTestType = "FISHER", correction = "BONFERRONI", type = "enrichment", panther_api.url = "http://pantherdb.org//services/oai/pantherdb/enrich/overrep", fdr.thresh = 0.1, p.thresh = 0.05){

  ## translate ontology and taxon into taxonID and ontology term
  oracl.settings = oracl_settings(ontology, taxon)

  ## build POST body
  oracl.body = oracl_body(x, bg, oracl.settings[1], oracl.settings[2], enrichmentTestType, correction, type)

  ## POST


  ## filter result

}
