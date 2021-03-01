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
#' @import dplyr
#' @param x A list of gene identifiers (character vector)
#' @param bg A list of background genes (character vector). If not provided, uses PANTHER standard background.
#' @param taxon Name of organism to be studied (Defaults to: Athaliana (ID: 3702))
#' @param ontology Panther annotation data type of interest (Defaults to "GO:0008150" for BP, also possible: "GO:0003674" for MF, "GO:0005575" for CC)
#' @param enrichmentTestType Type of enrichment test (Defaults to: "FISHER", also possible: BINOMIAL)
#' @param correction FP-Correction method (Defaults to: "FDR", also possible: BONFERRONI)
#' @param panther_api.url URL to API "http://pantherdb.org//services/oai/pantherdb/enrich/overrep"
#' @param fdr.thresh FDR cutoff (default: 0.1)
#' @param p.thresh P-value cutoff to be applied (default: p < 0.05)
#' @param only_overrepresented Retain only overrepresented GOterms? (default: TRUE)
#' @return A dataframe or panther output
#' @export

oraclient <- function(x,
                      bg = NULL,
                      ontology = "bp",
                      taxon = "Athaliana",
                      enrichmentTestType = "FISHER",
                      correction = "FDR",
                      panther_api.url = "http://pantherdb.org//services/oai/pantherdb/enrich/overrep",
                      fdr.thresh = 0.1,
                      p.thresh = 0.05,
                      only_overrepresented = TRUE){

  ## translate ontology and taxon into taxonID and ontology term
  oraclient.settings = oraclient_settings(ont = ontology, tax = taxon)
  fdr.thresh = ifelse(correction == "BONFERRONI", 0, fdr.thresh)

  ## build POST body
  oraclient.body = oraclient_body(geneset = x,
                          geneset_bg = bg,
                          ont = oraclient.settings[1],
                          tax = oraclient.settings[2],
                          ett = enrichmentTestType,
                          correct = correction)

  ## helpful message
  sprintf("[oraclient] Overrepresentation test. Ontology: %s, taxon: %s, %s test and %s correction. P-threshold is %s, FDR threshold is %s.",
          ontology,
          taxon,
          enrichmentTestType,
          correction,
          p.thresh,
          fdr.thresh) %>%
    print()

  ## POST
  oraclient.POST = oraclient_POST(oraclient.body, panther_api.url)
  oraclient.genes = oraclient.POST$results$input_list$mapped_id %>% do.call("c", .)

  ## get dataframe
  out.df = oraclient.POST %>%
    oraclient_to_df() %>%
    filter(fdr < fdr.thresh) %>%
    filter(p.value < p.thresh)

  ## add genes
  out_and_genes.df = oraclient_add_genes(df = out.df, genes = oraclient.genes, ont = ontology) %>%
    mutate_if(is.factor, as.character)

  ## format output
  if(only_overrepresented){

    out_and_genes.df = out_and_genes.df %>% filter(dir == "+")

  }

  return(out_and_genes.df)

}

#' Translate ontology and taxon to IDs
#'
#' translate ontology and taxon to IDs readable by the panther API.
#' @param ont User-specified ontology, one of "bp", "mf" or "cc".
#' @param tax User-specified taxon, must be "Athaliana".
#' @return A vector of length 2 with Panther-specific IDs for taxon and ontology.

oraclient_settings <- function(ont, tax){

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

oraclient_body <- function(geneset, geneset_bg = NULL, ont, tax, ett, correct){

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

oraclient_POST <- function(post_body, post_url){

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

#' Oraclient to dataframe
#'
#' Build a dataframe from output of oracl::oraclient_POST()
#' @param oraclient.ls A list of goterms
#' @return An ORA dataframe (tibble)
oraclient_to_df = function(oraclient.ls){

  goterms.ls = oraclient.ls$results$result
  goterms = lapply(goterms.ls, oraclient_extract_goterms) %>% do.call("rbind", .)

}

#' Extract GOterms from Panther output
#'
#' Extracts GOterms from Panther response
#' @param als Each and every entry in results$result of a Panther Enrichment
#' @return Tabular GOterms to be processed further
oraclient_extract_goterms = function(als){

  als.go = als$term
  goterm = ifelse(!is.null(als.go$id), yes = als.go$id, no = NA)
  label = ifelse(!is.null(als.go$label), yes = als.go$label, no = NA)

  data.frame("GO_ID" = goterm,
             "label" = label,
             "N" = als$number_in_list,
             "N.expected" = als$expected,
             "N.reference" = als$number_in_reference,
             "p.value" = als$pValue,
             "fdr" = als$fdr,
             "fold_change" = als$fold_enrichment,
             "dir" = als$plus_minus
  )
}

#' Add Gene Identifiers
#'
#' Add Locus IDs for overrepresented GOterms
#' @param als Each and every entry in results$result of a Panther Enrichment
#' @return Tabular GOterms to be processed further
oraclient_add_genes = function(df, genes, ont){

  if(nrow(df) != 0){
    ont.dict = list(oracl:::bp.df, oracl:::mf.df, oracl:::cc.df) %>%setNames(c("bp", "mf", "cc"))
    gos = df %>% pull(GO_ID) %>% as.character %>% na.omit
    ont.df = ont.dict[[ont]] %>% filter(locus %in% genes)

    ## get gene IDs per GOterm
    genes_ora = lapply(gos, grep, ont.df[[2]]) %>%
      setNames(gos) %>%
      lapply(function(x, y){

        out = y[x] %>% paste(collapse = ";")
        return(out)

      }, y = ont.df$locus)

    ## add gene ids to each GOterm
    gos_and_genes = lapply(names(genes_ora), function(x,y){

      term = x
      genelist = y[[x]]
      out = data.frame("GO_ID" = term,
                       "locus" = genelist)

    }, y = genes_ora) %>%
      do.call("rbind", .)

    ## add geneIds to GOterm-dataframe
    out = df %>% left_join(gos_and_genes)
  } else{

    out = df %>% mutate("locus" = "a")

  }

  ## output
  return(out)

}
