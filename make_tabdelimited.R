parseASJP <- function(fn, outfile) {
	if ( "stringr" %in% installed.packages()==FALSE ) {
		install.packages("stringr", repos="http://cran.us.r-project.org", verbose=FALSE, quiet=TRUE)
	}
	suppressMessages(library(stringr))
	if ( "MASS" %in% installed.packages()==FALSE ) {
		install.packages("MASS", repos="http://cran.us.r-project.org", verbose=FALSE, quiet=TRUE)
	}
	suppressMessages(library(MASS))
	db <- readLines(fn)
	#get rid of the string \tRosetta, which will upset things
	w_ros <- grep("\tRosetta", db)
	if ( length(w_ros) > 0 ) {
		for (i in 1:length(w_ros)) {
			db[w_ros[i]] <- strsplit(db[w_ros[i]], "\tRosetta")[[1]][1]
		}
	}
	#make a vector with the second lines of the metadata
	w_2nd <- grep("^ [[:digit:]] ", db)
	metadata2nd <- db[w_2nd]
	#make a vector with the first lines of the metadata
	w_1st <- w_2nd - 1
	metadata1st <- db[w_1st]
	#get rid of some exceptional cases of tabs in the first line of metadata
	if ( length(w_1st) > 0 ) {
		w_tab <- grep("\t", db[w_1st])
		if ( length(w_tab) > 0 ) {
			for (i in 1:length(w_tab)) {
				db[w_1st[w_tab[i]]] <- gsub("\t", "", db[w_1st[w_tab[i]]])
			}
		}
	}
	w_start <- w_1st[1]
	#in earlier versions of the database there could be multiple tabs
	#reduce multiple tabs to just one
	for (i in 1:length(db)) {
		gsub("\t.*", "\t", db[i])
	}
	# get rid of dots in datalines
	dataline_indices <- c(w_start:length(db))[-c(w_1st, w_1st+1)]
	for (i in 1:length(dataline_indices)) {
		db[dataline_indices[i]] <- gsub("\\.", "", db[dataline_indices[i]])
	}

	#extract metadata from metadata1st
	if ( length(metadata1st) > 0 ) {
		names <- c(); wls_fam <- c(); wls_gen <- c(); e <- c(); glot <- c()
		if ( (length(grep("\\{", metadata1st[1])) > 0 & length(grep("|", metadata1st[1])) > 0) & length(grep("@", metadata1st[1])) > 0 ) {
			for (i in 1:length(metadata1st)) {
				#extract names
				split <- strsplit(metadata1st[i],"\\{")
				names[i] <- split[[1]][1]
				#extract WALS family
				split <- strsplit(metadata1st[i],"\\{")
				split <- strsplit(split[[1]][2],"\\.")
				wls_fam[i] <- split[[1]][1]
				#extract WALS genus
				split <- strsplit(metadata1st[i],"\\.")
				split <- strsplit(split[[1]][2],"\\|")
				wls_gen[i] <- split[[1]][1]
				#extract Ethnologue classification
				split <- strsplit(metadata1st[i],"\\{")
				split <- strsplit(split[[1]][2],"\\|")
				split <- strsplit(split[[1]][2],"@")
				e[i] <- split[[1]][1]
				#extract Hammarström classification
				split <- strsplit(metadata1st[i],"\\{")
				split <- strsplit(split[[1]][2],"\\|")
				split <- strsplit(split[[1]][2],"@")
				split <- strsplit(split[[1]][2],"\\}")
				glot[i] <- split[[1]][1]
			}
		} else if ( length(grep("\\{", metadata1st[1])) > 0 & length(grep("|", metadata1st[1])) == 0 ) {
			for (i in 1:length(metadata1st)) {
				#extract names
				split <- strsplit(metadata1st[i],"\\{")
				names[i] <- split[[1]][1]
				#extract WALS family
				split <- strsplit(metadata1st[i],"\\{")
				split <- strsplit(split[[1]][2],"\\.")
				wls_fam[i] <- split[[1]][1]
				#extract WALS genus
				split <- strsplit(metadata1st[i],"\\.")
				split <- strsplit(split[[1]][2],"\\}")
				wls_gen[i] <- split[[1]][1]
			}
			e <- rep(NA, length(names))
			glot <- rep(NA, length(names))
		} else if ( (length(grep("\\{", metadata1st[1])) > 0 & length(grep("|", metadata1st[1])) > 0) & length(grep("@", metadata1st[1])) == 0 ) {
			for (i in 1:length(metadata1st)) {
				#extract names
				split <- strsplit(metadata1st[i],"\\{")
				names[i] <- split[[1]][1]
				#extract WALS family
				split <- strsplit(metadata1st[i],"\\{")
				split <- strsplit(split[[1]][2],"\\.")
				wls_fam[i] <- split[[1]][1]
				#extract WALS genus
				split <- strsplit(metadata1st[i],"\\.")
				split <- strsplit(split[[1]][2],"\\|")
				wls_gen[i] <- split[[1]][1]
				#extract Ethnologue classification
				split <- strsplit(metadata1st[i],"\\{")
				split <- strsplit(split[[1]][2],"\\|")
				split <- strsplit(split[[1]][2],"\\}")
				e[i] <- split[[1]][1]
			}
			glot <- rep(NA, length(names))
		} else {
			chunk <- db[w_start:length(db)]
			w_1stchunk <- grep("^ [[:digit:]] ", chunk) - 1
			names <- as.vector(unlist(lapply(chunk[w_1stchunk], function(s) strsplit(s, " ")[[1]][1])))
			wls_fam <- rep(NA, length(names))
			wls_gen <- rep(NA, length(names))
			e <- rep(NA, length(names))
			glot <- rep(NA, length(names))
			rm(chunk)
		}
	}
	#extract latitudes, longitudes, populations, WALS-codes, iso-codes from metadata2nd
	lat <- c(); lon <- c(); pop <- c(); wcode <- c(); iso <- c()
	for (i in 1:length(metadata2nd)) {
		lat[i] <- substr(metadata2nd[i],4,10)
		lon[i] <- substr(metadata2nd[i],12,18)
		pop[i] <- substr(metadata2nd[i],20,30)
		wcode[i] <- substr(metadata2nd[i],34,36)
		iso[i] <- substr(metadata2nd[i],40,42)
	}
	#delete all blanks
	names <- gsub(" ","",names)
	wls_fam <- gsub(" ","",wls_fam)
	wls_gen <- gsub(" ","",wls_gen)
	e <- gsub(" ","",e)
	glot <- gsub(" ","",glot)
	lat <- gsub(" ","",lat)
	lon <- gsub(" ","",lon)
	pop <- gsub(" ","",pop)
	wcode <- gsub(" ","",wcode)
	iso <- gsub(" ","",iso)
	#clean up
	rm(metadata1st)
	rm(metadata2nd)
	#at this point there are 10 vectors:
	#   names, wls_fam, wls_gen, e, glot, lat, lon, pop, wcode, iso
	#these are now combined into a dataframe
	dbf <- data.frame(names, wls_fam, wls_gen, e, glot, lat, lon, pop, wcode, iso, stringsAsFactors=FALSE)
	#clean up
	rm(names); rm(wls_fam); rm(wls_gen); rm(e); rm(glot); rm(lat); rm(lon); rm(pop); rm(wcode); rm(iso)
	#extract wordlists
	#fill columns with XXX beforehand
	v <- rep("XXX",length(dbf$names))
	dbf$I <- v; dbf$you <- v; dbf$we <- v; dbf$this <- v; dbf$that <- v; dbf$who <- v; dbf$what <- v; dbf$not <- v; dbf$all <- v; dbf$many <- v; dbf$one <- v; dbf$two <- v; dbf$big <- v; dbf$long <- v; dbf$small <- v; dbf$woman <- v; dbf$man <- v; dbf$person <- v; dbf$fish <- v; dbf$bird <- v; dbf$dog <- v; dbf$louse <- v; dbf$tree <- v; dbf$seed <- v; dbf$leaf <- v; dbf$root <- v; dbf$bark <- v; dbf$skin <- v; dbf$flesh <- v; dbf$blood <- v; dbf$bone <- v; dbf$grease <- v; dbf$egg <- v; dbf$horn <- v; dbf$tail <- v; dbf$feather <- v; dbf$hair <- v; dbf$head <- v; dbf$ear <- v; dbf$eye <- v; dbf$nose <- v; dbf$mouth <- v; dbf$tooth <- v; dbf$tongue <- v; dbf$claw <- v; dbf$foot <- v; dbf$knee <- v; dbf$hand <- v; dbf$belly <- v; dbf$neck <- v; dbf$breast <- v; dbf$heart <- v; dbf$liver <- v; dbf$drink <- v; dbf$eat <- v; dbf$bite <- v; dbf$see <- v; dbf$hear <- v; dbf$know <- v; dbf$sleep <- v; dbf$die <- v; dbf$kill <- v; dbf$swim <- v; dbf$fly <- v; dbf$walk <- v; dbf$come <- v; dbf$lie <- v; dbf$sit <- v; dbf$stand <- v; dbf$give <- v; dbf$say <- v; dbf$sun <- v; dbf$moon <- v; dbf$star <- v; dbf$water <- v; dbf$rain <- v; dbf$stone <- v; dbf$sand <- v; dbf$earth <- v; dbf$cloud <- v; dbf$smoke <- v; dbf$fire <- v; dbf$ash <- v; dbf$burn <- v; dbf$path <- v; dbf$mountain <- v; dbf$red <- v; dbf$green <- v; dbf$yellow <- v; dbf$white <- v; dbf$black <- v; dbf$night <- v; dbf$hot <- v; dbf$cold <- v; dbf$full <- v; dbf$new <- v; dbf$good <- v; dbf$round <- v; dbf$dry <- v; dbf$name <- v
	#identify all lines in the database where there are language names
	#then identify where there are words for the various concepts
	#w_n is a vector of where language names are
	#w_c is a list of vectors of where the concepts are
	#for each concept find the number in w_n which is closest
	#put the word for the concept in the matrix
	#at the row corresponding to the language name and the 
	#column corresponding to the concept
	w_n <- grep("^ [[:digit:]] ", db) - 1
	w_c <- vector(mode = "list")
	for (i in 1:100) {
		w_c[[i]] <- c(grep(paste("^0*",i," ",sep=""),db))
	}
	for (i in 1:100) {
#		cat("Processing Swadesh list item ", i,"\n")
		for (j in 1:length(w_c[[i]])) {
			ord <- sort(c(w_n,w_c[[i]]))
			lg_no <- which(ord[which(w_c[[i]][j]==ord)-1]==w_n)
			item <- str_extract(db[w_c[[i]][j]],"\t.*/*")
			item <- gsub(" /","/",item)
			item <- gsub("/.*","",item)
			item <- gsub("\t","",item)
			dbf[lg_no,i+10] <- item
		}
	}
	write.table(dbf,file=outfile,row.names=FALSE,quote=FALSE,sep="\t")
}
