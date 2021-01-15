###############################
#
# Generate GOterm datasets
#
###############################

bp.df = read.csv("data/goterms/go_bp_ath.csv", stringsAsFactors = FALSE)
mf.df = read.csv("data/goterms/go_mf_ath.csv", stringsAsFactors = FALSE)
cc.df = read.csv("data/goterms/go_cc_ath.csv", stringsAsFactors = FALSE)

usethis::use_data(bp.df, mf.df, cc.df, internal = TRUE)
