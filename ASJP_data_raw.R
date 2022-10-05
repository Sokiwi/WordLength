# load required functions
source("word_length_entry.R")  # loads wl() function
source("synonyms_entry.R")  # loads syn() function

# put the current version of the ASJP database
# in a tab-delimited format, producing listss19_formatted.tab
# or listss20_formatted.tab etc., depending on the version number;
# current version is called lists.txt and is found in the zipped file
# at https://zenodo.org/record/3840975
# later versions will be in an equivalent place
# only needed if listss19_formatted.tab (or similar) is not available
# load parseASJP() function:
source("make_tabdelimited.R")
# run the formatter, indicating in- and output names
# this takes a while, but less than an hour
parseASJP("listss19.txt", "listss19_formatted.tab")

## Prepare ASJP data raw
# modify the name of the input file as needed
x <- read.table(file="listss19_formatted.tab", sep="\t", 
   quote="",stringsAsFactors=FALSE, header=TRUE, na.strings="")

# exclude doculects that have PROTO in their name or are tagged 
# as ARTIFICIAL, FAKE or SPEECH_REGISTER in their WALS 
# genus classification
excl1 <- grep("PROTO", x$names)
excl2 <- which(x$wls_gen=="ARTIFICIAL")
excl3 <- which(x$wls_gen=="FAKE")
excl4 <- which(x$wls_gen=="SPEECH_REGISTER")
# exclude doculects with no Glottolog classification in the dataset
excl5 <- which(x$glot=="NA")
# exclude doculects classified as MixedLanguage, Spurious, Unclassifiable
# or Artificial in Glottolog
excl6 <- which(x$glot %in% c("MixedLanguage", "Spurious", 
  "Unclassifiable_xil", "Unclassifiable_kzw", "ArtificialLanguage"))
# exclude doculects only documented by the year 1700 or earlier
excl7 <- which(x$pop >= -1700 & x$pop <= 0)
excl <- unique(c(excl1, excl2, excl3, excl4, excl5, excl6, excl7))
x <- x[-excl,]

# make a matrix and fill it with word length counts
# (takes around 1 min)
xwl <- matrix(NA, nrow=nrow(x), ncol=100)
total <- nrow(x)*100
count <- 0
for (i in 1:nrow(x)) {
	for (j in 11:110) {
		count <- count + 1
		xwl[i,j-10] <- wl(x[i,j])
		if ( count%%1000 == 0 ) {
			cat("working on", count, "out of", total, "\n")
		}
	}
}

# make a matrix and fill it with counts of number of synonyms
# (takes around 1 min)
xsyn <- matrix(NA, nrow=nrow(x), ncol=100)
total <- nrow(x)*100
count <- 0
for (i in 1:nrow(x)) {
	for (j in 11:110) {
		count <- count + 1
		xsyn[i,j-10] <- syn(x[i,j])
		if ( count%%1000 == 0 ) {
			cat("working on", count, "out of", total, "\n")
		}
	}
}

# aggregate by ISO codes
# If there is no ISO code the ASJP name will be pasted together with
# the string "no_iso_" and will serve as key when aggregating;
# If the ISO code is a code supplied by ASJP for a language that
# really should have an ISO code it will be pasted together with
# "no_iso_". There is least one case of more than one list with the
# same pseudo-ISO code
for (i in 1:nrow(x)) {
	if ( is.na(x$iso[i]) ) {
		x$iso[i] <- paste("no_iso_", x$names[i], sep="")
	}
	if ( length(grep("[[:lower:]][[:lower:]]0", x$iso[i])) == 1 ) {
		x$iso[i] <- paste("no_iso_", x$iso[i], sep="")
	}
}

# make a data frame with the meta data, as outlined in the data
# preparation document. The names column will contain the ASJP names,
# all other metadata will be selected from the first ISO-code
# member
iso <- unique(x$iso)
names <- c()
for (i in 1:length(iso)) {
	w <- which(x$iso==iso[i])
	names[i] <- paste(x$names[w], collapse=",")
}
# identify first occurrence of each unique iso-code
# for the purpose of extracting meta data
fo <- match(iso, x$iso)  # first occurrence
# prepare columns of meta data for WALS code, WALS family,
# Glottolog fammily, WALS genus, latitude, longitude, population
wals_code <- x$wcode[fo]
wals_fam <- x$wls_fam[fo]
get_fam <- function(z) {
	strsplit(z, ",")[[1]][1]
}
glot_fam_all <- unlist(lapply(x$glot, get_fam))
glot_fam <- glot_fam_all[fo]
wals_genus <- x$wls_gen[fo]
lat <- x$lat[fo]
lon <- x$lon[fo]
pop <- x$pop[fo]
for (i in 1:length(pop)) {
	if ( pop[i] < 1 ) {
		pop[i] <- NA
	}
}

# aggregate wl values for ISO-codes
xwl_agg <- matrix(NA, nrow=length(iso), ncol=100)
for (i in 1:length(iso)) {
	w_iso <- which(x$iso==iso[i])
	for (j in 1:100) {
		wl_agg <- round(mean(xwl[w_iso,j], na.rm=TRUE), 3)
		if (is.na(wl_agg)) {
			wl_agg <- NA
		}
		xwl_agg[i,j] <- wl_agg
	}
}

# aggregate synonym counts for ISO-codes
xsyn_agg <- matrix(NA, nrow=length(iso), ncol=100)
for (i in 1:length(iso)) {
	w_iso <- which(x$iso==iso[i])
	for (j in 1:100) {
		syn_agg <- round(100*(mean(xsyn[w_iso,j], na.rm=TRUE) - 1), 2)
		if (is.na(syn_agg)) {
			syn_agg <- NA
		}
		xsyn_agg[i,j] <- syn_agg
	}
}

# turn xwl_agg into a data frame with column, to be merged with
# the meta data
xwl_agg <- as.data.frame(xwl_agg)
names(xwl_agg) <- names(x)[11:110]

# turn xsyn_agg into a data frame with column, to be merged with
# the meta data
xsyn_agg <- as.data.frame(xsyn_agg)
names(xsyn_agg) <- paste("syn_", names(x)[11:110], sep="")

# prepare a column just containing either 40 or 100 in order to 
# signal the nature of the wordlist; here 100 would be used if 
# at least one of the lists aggregated over is a 100-item list
# (this will precede the previous data column with word counts in
# the final data frame)

# define the the non-40 and the 40 items
# the non-forty were defined earlier in thie script, but on a data frame 
# with some column excluded, so this definition is taken but moved by 9
non_forty <- c(5:11, 14:18, 21, 25, 27:28, 30, 33:34, 36:39, 43, 46:47, 
  50:51, 53, 56:57, 60:61, 63:66, 68:72, 74, 77, 79:82, 84:85, 88:92, 
  94:95, 98:100)
non_forty <- non_forty - 1
forty <- setdiff(1:100, non_forty)
forty_hundred <- c()
for (i in 1:nrow(xwl_agg)) {
	if ( sum(xwl_agg[i,non_forty], na.rm=TRUE)==0 ) {
		forty_hundred[i] <- "F"
	} else {
		forty_hundred[i] <- "H"
	}
}

# prepare column for number of attestations
perc_att_40 <- c()
perc_att_100 <- c()
for (i in 1:nrow(xwl_agg)) {
	missing_40 <- length(which(is.na(xwl_agg[i,forty])))
	perc_att_40[i] <- round(100*(40 - missing_40)/40,2)
	missing_100 <- length(which(is.na(xwl_agg[i,])))
	perc_att_100[i] <- 100 - missing_100
}

# get mean word length over the 40 and over all 100 items
forty_mean <- c()
hundred_mean <- c()
for (i in 1:nrow(xwl_agg)) {
	forty_mean_value <- round(mean(as.numeric(xwl_agg[i,forty]), na.rm=TRUE),3)
	if (is.na(forty_mean_value)) {
		forty_mean_value <- NA
	}
	forty_mean[i] <- forty_mean_value
	hundred_mean_value <- round(mean(as.numeric(xwl_agg[i,]), na.rm=TRUE),3)
	if (is.na(hundred_mean_value)) {
		hundred_mean_value <- NA
	}
	hundred_mean[i] <- hundred_mean_value
}
nas <- which(forty_hundred=="F")
hundred_mean[nas] <- NA
ASJP_data_raw <- data.frame(iso, names, wals_code, wals_fam, glot_fam, 
  wals_genus, lat, lon, pop, forty_hundred, xwl_agg, xsyn_agg, perc_att_40,
  perc_att_100, forty_mean, hundred_mean)

write.table(ASJP_data_raw, file="Data-01 ASJP data raw.txt", sep="\t", 
  quote=FALSE, row.names=FALSE)

# add AutoTyp areas to ASJP_data_raw
# use the line below in case ASJP_data_raw is not in  memory
ASJP_data_raw <- read.table(file="Data-01 ASJP data raw.txt", header=TRUE, sep="\t", quote="")

# before reading the Autotyp data it was reduced to relevant columns
# languages associated with errors as per file errors_Autotyp.txt were removed
# also Koyra had wrong assignment and was removed
# languages without ISO-codes were removed
# languages without area information was removed
# languages without coordinates were removed
# 2898 (Cajun French) was removed because the coordinates would misplace French (fra)
# 1406 (Spanish Canary Islands) was removed because the coordinates would misplace Spanish (spa)
# 23 (Arabic) was removed because the coordinates would put Classic Arabic in North Africa
# a pseudo-member of Autotyp was added (LID 100000, jbn), given correct metadata
# as a trick to avoid that it ends up in Europe;
# a new file formatted like the one supplied should be prepared for later versions of AutoTyp
AT <- read.csv(file="autotyp_area_data.csv", stringsAsFactors=FALSE, comment.char="", quote="", 
  header=TRUE, na.strings="NA")  # the AutoType data

unique_isos <- unique(AT$ISO639.3)
w_u_i <- match(unique_isos, AT$ISO639.3)
AT <- AT[w_u_i,]

# add areas and macroareas for ISO-codes that are
# also in AutoTyp (AT); if they aren't take the data
# for the closest language
# (takes around a minute)
source("distance_function.R")  # contains the distance() function from defunct argosfilter
area <- c()
continent <- c()
macrocontinent <- c()
distances <- c()
nearest_language <- c()
for (i in 1:nrow(ASJP_data_raw)) {
	w_iso <- which(AT$ISO639.3==ASJP_data_raw$iso[i])
	if ( length(w_iso) > 0 ) {
		area[i] <- AT$Area[w_iso]
		continent[i] <- AT$Continent[w_iso]
		macrocontinent[i] <- AT$Macrocontinent[w_iso]
		distances[i] <- 0
		nearest_language[i] <- "iso_twin"	
	} else {
		d <- c()  # a vector of distances to all languages
		for (j in 1:nrow(AT)) {
			d[j] <- distance(ASJP_data_raw$lat[i], AT$Latitude[j], 
                     ASJP_data_raw$lon[i], AT$Longitude[j])
		}
		w_min <- which(d==min(d))[1]
		area[i] <- AT$Area[w_min]
		continent[i] <- AT$Continent[w_min]
		macrocontinent[i] <- AT$Macrocontinent[w_min]
		distances[i] <- d[w_min]
		nearest_language[i] <- AT$ISO639.3[w_min]	
	}
}

ASJP_data_raw <- data.frame(ASJP_data_raw, area, continent, macrocontinent)

write.table(ASJP_data_raw, file="Data-01 ASJP data raw.txt", sep="\t", 
  quote=FALSE, row.names=FALSE)

