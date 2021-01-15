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
#' @param bg A list of background genes (character vector). If not provided, uses PANTHER standard background.
#' @param taxon Name of organism to be studied (Defaults to: Athaliana (ID: 3702))
#' @param ontology Panther annotation data type of interest (Defaults to "GO:0008150" for BP, also possible: "GO:0003674" for MF, "GO:0005575" for CC)
#' @param enrichmentTestType Type of enrichment test (Defaults to: "FISHER", also possible: BINOMIAL)
#' @param correction FP-Correction method (Defaults to: "BONFERRONI", also possible: FDR)
#' @param panther_api.url URL to API "http://pantherdb.org//services/oai/pantherdb/enrich/overrep"
#' @param fdr.thresh FDR cutoff (default: 0.1)
#' @param p.thresh P-value cutoff to be applied (default: p < 0.05)
#' @return A dataframe
#' @export

oracl <- function(x, bg = NULL, ontology = "bp", taxon = "Athaliana",  enrichmentTestType = "FISHER", correction = "BONFERRONI", panther_api.url = "http://pantherdb.org//services/oai/pantherdb/enrich/overrep", fdr.thresh = 0.1, p.thresh = 0.05){

  ## translate ontology and taxon into taxonID and ontology term
  oracl.settings = oracl_settings(ont = ontology, tax = taxon)

  ## build POST body
  oracl.body = oracl_body(geneset = x,
                          geneset_bg = bg,
                          ont = oracl.settings[1],
                          tax = oracl.settings[2],
                          ett = enrichmentTestType,
                          correct = correction)

  ## POST
  oracl.POST = oracl_POST(oracl.body, panther_api.url)

  ## filter result

}

#' Translate ontology and taxon to IDs
#'
#' translate ontology and taxon to IDs readable by the panther API.
#' @param ont User-specified ontology, one of "bp", "mf" or "cc".
#' @param tax User-specified taxon, must be "Athaliana".
#' @return A vector of length 2 with Panther-specific IDs for taxon and ontology.

oracl_settings <- function(ont, tax){

  ## prepate dictionaries
  ont.dict = c("GO:0008150", "GO:0003674", "GO:0005575") %>%  setNames(c("bp", "mf", "cc"))
  tax.dict = 3702 %>% setNames("Athaliana")

  ## output mapped IDs
  return(
    c(ont.dict[[ont]],
      tax.dict[[tax]])
  )

}

#' Translate ontology and taxon to IDs
#'
#' translate ontology and taxon to IDs readable by the panther API.
#' @param geneset User-specified geneset
#' @param geneset_bg Background. If null, use Panther standard
#' @param ont ontology to be investigated (ID!)
#' @param tax Taxon (ID!)
#' @param ett enrichmentTestType
#' @param correct correction method
#' @param tax User-specified taxon, must be "Athaliana".
#' @return A vector of length 2 with Panther-specific IDs for taxon and ontology.

oracl_body <- function(geneset, geneset_bg = NULL, ont, tax, ett, correct){

  ## if background is specified, build body with background geneset
  if(!is.null(geneset_bg)){

    out = list(
      "organism" = tax,
      "refOrganism" = tax,
      "refInputList" = geneset_bg %>% paste(collapse = ","),
      "geneInputList" = geneset %>% paste(collapse = ","),
      "annotDataSet" = ont,
      "enrichmentTestType" = ett,
      "correction" = correct
    )

  }

  ## else build without background geneset
  else {

      out = list(
        "organism" = tax,
        "geneInputList" = geneset %>% paste(collapse = ","),
        "annotDataSet" = ont,
        "enrichmentTestType" = ett,
        "correction" = correct
      )
  }

  ## get result
  return(out)

}

#' POST request
#'
#' Use predefined body and API-URL for post request
#' @importFrom httr POST
#' @importFrom httr content
#' @importFrom retry retry
#' @param post_body Body to be used (a list)
#' @param post_url URL to be used (character)
#' @return An ORA (.json)

oracl_POST <- function(post_body, post_url){

  ## POST request and extraction of content
  out = retry(

    POST(

      url = post_url,
      body = post_body,
      encode = "form"

    ),
    when = "Connection was reset"

  ) %>%
    content()

  ## spuck aus
  return(out)

}

