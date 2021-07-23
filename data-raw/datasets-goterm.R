###############################
#
# Generate GOterm datasets
#
###############################

## GO-Term datasets
bp.df = read.csv("data/goterms/go_bp_ath.csv", stringsAsFactors = FALSE)
mf.df = read.csv("data/goterms/go_mf_ath.csv", stringsAsFactors = FALSE)
cc.df = read.csv("data/goterms/go_cc_ath.csv", stringsAsFactors = FALSE)

## sample datasets
GS01 = readLines("data/genesets/GS01.txt")
GS02 = readLines("data/genesets/GS02.txt")
GS03 = readLines("data/genesets/GS03.txt")
GS04 = readLines("data/genesets/GS04.txt")
GS05 = readLines("data/genesets/GS05.txt")
background = readLines("data/background/background.txt")

## export
usethis::use_data(GS01,
                  GS02,
                  GS03,
                  GS04,
                  GS05,
                  background,
                  bp.df,
                  mf.df,
                  cc.df,
                  internal = TRUE,
                  overwrite = TRUE)
