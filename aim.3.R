library(plyr) # load before dplyr to avoid conflicts due to masking
library(dplyr)
library(tcR)
library(ggplot2)
library(grid)
library(gridBase)
library(ggrepel)
library(magrittr)
library(stringi)
library(stringr)
library(DescTools)
source("tcR_functions.R") #includes several modifications to tcR

switch(Sys.info()[['user']],
       nealp = {fig.file.path <- "C:/Users/nealp/Dropbox (Personal)/Extension School/Thesis/figures"
       raw.file.path <- "C:/Users/nealp/Dropbox (Partners HealthCare)/Projects/PNOIT2-1037/TCRB sequencing and HLA typing data/neals.thesis.data"
       subset.data.path = "C:/Users/nealp/Dropbox (Partners HealthCare)/PAID-UP, PNOIT2 Reactive vs Non-reactive paper/TCR Sequencing/Raw data BL T cell subsets"},
       wayne1 = {fig.file.path <- "~/Dropbox (Personal)/Extension School/Thesis/figures"
       raw.file.path <- "~/Dropbox (Partners HealthCare)/Projects/PNOIT2-1037/TCRB sequencing and HLA typing data/neals.thesis.data"
       subset.data.path = "~/Dropbox (Partners HealthCare)/PAID-UP, PNOIT2 Reactive vs Non-reactive paper/TCR Sequencing/Raw data BL T cell subsets"},
       ns580 = {fig.file.path <- "~/thesis/figures"
       raw.file.path <- "~/data/thesis"},
       stop("I don't recognize your username, type Sys.info() to find out what it is.")
)

switch(Sys.info()[['user']],
       nealp = {neals.func.path = "C:/Users/nealp/Dropbox (Partners HealthCare)/Projects/PNOIT2-1037/TCRB sequencing and HLA typing data/neals_tcr_functions.R"},
       wayne1 = {neals.func.path = "C:/Users/nealp/Dropbox (Partners HealthCare)/Projects/PNOIT2-1037/TCRB sequencing and HLA typing data/neals_tcr_functions.R"},
       ns580 = {neals.func.path = "~/data/neals_tcr_functions.R"},
       stop("I don't recognize your username, type Sys.info() to find out what it is.")
)
source(neals.func.path)


# Load in the data
load(paste(raw.file.path, "parsed.data.baseline.rda", sep = "/"))
load(paste(raw.file.path, "enriched.CDR3s.rda", sep = "/"))

enriched.CDR3s.df <- do.call(rbind, enriched.CDR3s)


# Includes 4 hyporeactives who were confirmed allergic by challenge!
reactive.ids <- c("69", "80", "84", "85", "89", "93", "100", "104", "106", "107", "27",
                  "81", "101", "94")

# Get the files
reactive.subset.files <- list.files(paste(sub("/[^/]*$", "", raw.file.path), "teff.treg.data", sep = "/"))
hypo.conf.files <- list.files(paste(subset.data.path, "Hyporeactive", sep = "/"))[
  grep(paste(reactive.ids, sep = "", collapse = "|"), 
       list.files(paste(subset.data.path, "Hyporeactive", sep = "/")))]

# Parse all of the data
subset.parse <- parse.adaptive.file.list(c(paste(sub("/[^/]*$", "", raw.file.path), "/teff.treg.data", reactive.subset.files, sep = "/"),
                                           paste(subset.data.path, "/Hyporeactive", hypo.conf.files, sep = "/")))

# Seperate the data into different subsets
teff.parse <- subset.parse[grep("Teff", names(subset.parse))]
treg.parse <- subset.parse[grep("Treg", names(subset.parse))]

# Limit each to just the enriched CDR3s
teff.psCDR3 <- lapply(teff.parse, function(x){
  x[x$CDR3.amino.acid.sequence %in% enriched.CDR3s.df$CDR3.amino.acid.sequence,]
})

# Limit each to just the enriched CDR3s
treg.psCDR3 <- lapply(treg.parse, function(x){
  x[x$CDR3.amino.acid.sequence %in% enriched.CDR3s.df$CDR3.amino.acid.sequence,]
})

# Load in information about clusters
load(paste(raw.file.path, "graph.layout.rda", sep = "/"))

# Write a function to determine which clusters the CDR3s in the subset data are in
x <- lapply(teff.psCDR3, function(x){
  ids = c()
  for(i in 1:nrow(x)){
    if(x[i,]$CDR3.amino.acid.sequence %in% graph.layout$CDR3aa){
      clust <- graph.layout$cid[graph.layout$CDR3aa == x[i,]$CDR3.amino.acid.sequence]
      ids = c(clust, ids)
    }
  }
  return(ids)
})

