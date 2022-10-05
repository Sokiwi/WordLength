# function for counting the number of synonyms in an ASJP entry
# counts phrases same as single words, and returns 2 if there
# are 2 or more synonyms
syn <- function(z) {
	z1 <- gsub("XXX", "", z)
	if ( is.na(z) == TRUE ) {
		return(NA)
	} else if ( nchar(z1) == 0 ) {
		return(NA)
	} else {
	# make a list of synonyms
		syn1 <- strsplit(z, ",")[[1]]
		synonyms <- length(syn1)
		if ( synonyms > 2 ) {
			synonyms <- 2
		}
		return(synonyms)
	}
}
