# Resampling approach to determine number of edges in control data

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
library(stringdist)
library(RVAideMemoire)
library(DescTools)

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
load(paste(raw.file.path, "top.nmers.rda", sep = "/"))
load(paste(raw.file.path, "top.disc.threemers.rda", sep = "/"))

enriched.CDR3s.df <- do.call(rbind, enriched.CDR3s)

# Put top nmers in a dataframe
top.nmer.df <- do.call(rbind, top.nmers)
top.nmer.df$nmer <- as.character(top.nmer.df$nmer)
# For accurate neighbor count, need to make sure we get unique nucleotide sequences, not just AA sequences
pos.parse <- data.parse[grep("pos", names(data.parse))]
names(pos.parse) <- substr(names(pos.parse), 0, nchar(names(pos.parse)) - 4)

# # TEST: Make sure names are correct and in order
# names(pos.parse) == names(enriched.CDR3s)

for(i in 1:length(pos.parse)){
  pos.parse[[i]] <- pos.parse[[i]][pos.parse[[i]]$CDR3.amino.acid.sequence %in%
                                     enriched.CDR3s[[i]]$CDR3.amino.acid.sequence,]
}

# Get just the negative data
neg.parse <- data.parse[grep("neg", names(data.parse))]
names(neg.parse) <- substr(names(neg.parse), 0, nchar(names(neg.parse)) - 4)

# # TEST: Make sure this is true to ensure proper sampling
# names(neg.parse) == names(pos.parse)

# Function to find CDR3s that are within 1 AA difference (lev. distance)
find_pairs_hom <- function(x, y) {
  res <- dist.to.df(as.dist(stringdistmatrix(x, y,
                                             method = "lv",
                                             useNames = "strings"))) %>%
    filter(dist == 1) %>%
    select(-dist)
  colnames(res) <- c("from.cdr3", "to.cdr3")
  res
}


# Get the number of CDR3s to resample from each individual
num.CDR3s <- unlist(lapply(pos.parse, nrow), use.names = FALSE)

# Get the negative data
CD154.neg.CDR3s <- lapply(data.parse[grep("neg", names(data.parse))], function(x){
  x = select(x, CDR3.amino.acid.sequence, Read.count)
  return(x)
}) %>% do.call(rbind, .)

# Aggregate all repeated CDR3s, adding up read counts
CD154.neg.CDR3s <- aggregate(. ~ CDR3.amino.acid.sequence, data = CD154.neg.CDR3s, sum)

# Get all possible nmer data from the negative CDR3s: will make resampling faster and less redundent
# Create list of 3mers, 4mers and 5mers
threemer.list.neg <- nmer.gen(CD154.neg.CDR3s$CDR3.amino.acid.sequence, 3)
fourmer.list.neg <- nmer.gen(CD154.neg.CDR3s$CDR3.amino.acid.sequence, 4)
fivemer.list.neg <- nmer.gen(CD154.neg.CDR3s$CDR3.amino.acid.sequence, 5)

# Set cores for parallel processing
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

# Use paralell processing to look for the nmers in the CD154- CDR3s
threemer.neg <- rbindlist(foreach(i = 1:length(threemer.list.neg)) %dopar% {
  source(paste(raw.file.path, "neals_tcr_functions.R", sep = "/"))
  data = nmer.search(threemer.list.neg[i], CD154.neg.CDR3s, CDR3.col = 1, count.col = 2)
  data
})

# Set cores for parallel processing
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

# Use paralell processing to look for the nmers in the CD154- CDR3s
fourmer.neg <- rbindlist(foreach(i = 1:length(fourmer.list.neg)) %dopar% {
  source(paste(raw.file.path, "neals_tcr_functions.R", sep = "/"))
  data = nmer.search(fourmer.list.neg[i], CD154.neg.CDR3s, CDR3.col = 1, count.col = 2)
  data
})

# Set cores for parallel processing
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

# Use paralell processing to look for the nmers in the CD154- CDR3s
fivemer.neg <- rbindlist(foreach(i = 1:length(fivemer.list.neg)) %dopar% {
  source(paste(raw.file.path, "neals_tcr_functions.R", sep = "/"))
  data = nmer.search(fivemer.list.neg[i], CD154.neg.CDR3s, CDR3.col = 1, count.col = 2)
  data
})

# Put all of the data into a list
motif.neg <- list(threemer = threemer.neg,
                  fourmer = fourmer.neg,
                  fivemer = fivemer.neg)

# Save this list as an rda


# Set up resampling
for(resamp in 1:1){
  # get random CDR3s
  random.CDR3s <- c()
  random.CDR3.df <- list()
  # Get same amount of CDR3s from each individual that they contribute to our enricheed set
  for(i in 1:length(pos.parse)){
    CDR3s <- sample(neg.parse[[i]]$CDR3.amino.acid.sequence, num.CDR3s[i])
    random.CDR3.df[[i]] <- neg.parse[[i]][neg.parse[[i]]$CDR3.amino.acid.sequence %in% CDR3s,]
    random.CDR3s <- c(random.CDR3s, CDR3s)
  }
  #pairs.lev <- find_pairs_hom(random.CDR3s, random.CDR3s)
  
  # Motif analysis on random CDR3s
  random.CDR3.df <- do.call(rbind, random.CDR3.df)
  
  # Create list of 3mers, 4mers and 5mers
  threemer_list <- nmer.gen(random.CDR3.df$CDR3.amino.acid.sequence, 3)
  fourmer_list <- nmer.gen(random.CDR3.df$CDR3.amino.acid.sequence, 4)
  fivemer_list <- nmer.gen(random.CDR3.df$CDR3.amino.acid.sequence, 5)
  
  # Determine proportions of nmers in random CDR3s
  threemer.random <- nmer.search(threemer_list, random.CDR3.df, count.col = 3)
  fourmer.random <- nmer.search(fourmer_list, enriched.CDR3s.df, count.col = 3)
  fivemer.random <- nmer.search(fivemer_list, enriched.CDR3s.df, count.col = 3)
  
  # Set cores for parallel processing
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1)
  registerDoParallel(cl)
  
  # Use paralell processing to look for the nmers in the CD154- CDR3s
  threemer.neg <- rbindlist(foreach(i = 1:length(threemer.random$nmer)) %dopar% {
    source(paste(raw.file.path, "neals_tcr_functions.R", sep = "/"))
    data = nmer.search(threemer.random$nmer[i], CD154.neg.CDR3s, CDR3.col = 1, count.col = 2)
    data
  })
  
  # Set cores for parallel processing
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1)
  registerDoParallel(cl)
  
  fourmer.neg <- rbindlist(foreach(i = 1:length(fourmer.random$nmer)) %dopar% {
    source(paste(raw.file.path, "neals_tcr_functions.R", sep = "/"))
    data = nmer.search(fourmer.enriched$nmer[i], CD154.neg.CDR3s, CDR3.col = 1, count.col = 2)
    data
  })
  
  
}





