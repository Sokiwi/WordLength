# function for counting the word length for a an ASJP entry
wl <- function(z) {
	z1 <- gsub("XXX", "", z)
	if ( is.na(z) == TRUE ) {
		return(NA)
	} else if ( nchar(z1) == 0 ) {
		return(NA)
	# return NA when there is just one synonym and this is a phrase
	} else if ( length(grep(",", z))==0 & length(grep(" ", z)) > 0 ) {
		return(NA)
	} else {
	# make a list of synonyms
		syn1 <- strsplit(z, ",")[[1]]
		if ( length(syn1) > 2 ) {
			syn1 <- syn1[1:2]
		}
		# get rid of trailing spaces
		syn2 <- as.vector(unlist(lapply(syn1, trimws)))
		# assign NA values to phrases
		for (i in 1:length(syn2)) {
			if ( length(grep(" ", syn2[i])) > 0 ) {
				syn2[i] <- NA
			}
		}
		# return NA if there were only (several) phrases
		if ( length(which(is.na(syn2)))==length(syn2) ) {
			return(NA)
		} else {
			# count word length for each synonym
			for (j in 1:length(syn2)) {
				if ( !is.na(syn2[j]) ) {
					w <- syn2[j]
					w1 <- gsub("%", "", w)
					w2 <- gsub("\\*", "", w1)
					w3 <- gsub("[[:alnum:]]{2}~", "Q", w2)
					w4 <- gsub("[[:alnum:]]{3}\\$", "Q", w3)
					w5 <- gsub("\\\"", "", w4)
					syn2[j] <- w5
				}
			}
		}
		word_length <- round(mean(nchar(syn2), na.rm=TRUE),2)
		if ( word_length==0 ) {
			return(NA)
		} else {
			return(word_length)
		}
	}
}

