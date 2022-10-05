## make table called WALS_data.txt
# load required functions
# next line only needed if "listss19_formatted.tab" is not available
# source("make_tabdelimited.R")  # loads parseASJP() function
source("word_length_entry.R")  # loads wl() function

# run the formatter (only needed if "listss19_formatted.tab" 
# is not available)
# parseASJP("listss19.txt", "listss19_formatted.tab")
x <- read.table(file="listss19_formatted.tab", sep="\t", 
   quote="",stringsAsFactors=FALSE, header=TRUE, na.strings="")
# select only those that have WALS codes
x <- x[-which(is.na(x$wcode)),]
# remove unwanted columns, both metadata and entries not
# pertaining to the 40-item list
x <- x[,-c(1:8,10)]
non_forty <- c(5:11, 14:18, 21, 25, 27:28, 30, 33:34, 36:39, 43, 
  46:47, 50:51, 53, 56:57, 60:61, 63:66, 68:72, 74, 77, 79:82, 
  84:85, 88:92, 94:95, 98:100)
x <- x[,-non_forty]
# make a matrix and fill it with word length counts
# (takes around 1 min)
xwl <- matrix(NA, nrow=nrow(x), ncol=40)
total <- nrow(x)*(ncol(x)-1)
count <- 0
for (i in 1:nrow(x)) {
	for (j in 2:ncol(x)) {
		count <- count + 1
		xwl[i,j-1] <- wl(x[i,j])
		if ( count%%1000 == 0 ) {
			cat("working on", count, "out of", total, "\n")
		}
	}
}

# aggregate by WALS codes
codes <- sort(unique(x$wcode))
xwl_agg <- matrix(NA, nrow=length(codes), ncol=40)
row.names(xwl_agg) <- codes
for (i in 1:length(codes)) {
	wcode <- which(x$wcode==codes[i])
	for (j in 1:40) {
		wl_agg <- mean(xwl[wcode,j], na.rm=TRUE)
		if (is.na(wl_agg)) {
			wl_agg <- NA
		}
		xwl_agg[i,j] <- wl_agg
	}
}
means <- c()
for (i in 1:nrow(xwl_agg)) {
	means[i] <- round(mean(as.numeric(xwl_agg[i,1:40]), na.rm=TRUE),3)
}

# produce WALS data table with word length as a continuous
# variable, including only sufficiently attested languages 
# (20 or less missing items)
# and also including categorical values

nas <- apply(xwl_agg, 1, function(z) length(which(is.na(z))))
w_suf <- which(nas <= 20)  # languages with sufficient attestation
WALS_data <- data.frame(codes[w_suf],means[w_suf])
colnames(WALS_data) <- c("code","MWL")
value <- cut(WALS_data$MWL, breaks=c(0, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 8), 
  labels=c("extremely short", "very short", 
  "short", "moderately short", "moderately long", "long", "very long",
  "extremely long"))
WALS_data <- cbind(WALS_data, value)
write.table(WALS_data, file="Data-02 WALS data.txt", sep="\t", 
  quote=FALSE, row.names=FALSE)
