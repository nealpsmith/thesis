### TCRNet analysis ###
rm(list=ls())
library(plyr) # load before dplyr to avoid conflicts due to masking
library(dplyr)
library(tcR)
library(ggplot2)
library(GGally)
library(grid)
library(gridBase)
library(ggrepel)
library(magrittr)
library(stringdist)
library(network)
library(doParallel)
library(foreach)
library(DescTools)

source("tcR_functions.R") #includes several modifications to tcR
source("neals_tcr_functions.R")

switch(Sys.info()[['user']],
       nealp = {fig.file.path <- "C:/Users/nealp/Dropbox (Personal)/Extension School/Thesis/figures"
       raw.file.path <- "C:/Users/nealp/Dropbox (Partners HealthCare)/Projects/PNOIT2-1037/TCRB sequencing and HLA typing data/neals.thesis.data"
       subset.data.path = "C:/Users/nealp/Dropbox (Partners HealthCare)/PAID-UP, PNOIT2 Reactive vs Non-reactive paper/TCR Sequencing/Raw data BL T cell subsets"},
       wayne1 = {fig.file.path <- "~/Desktop/Neal_temp/figures"
       raw.file.path <- "~/Dropbox (Partners HealthCare)/Projects/PNOIT2-1037/TCRB sequencing and HLA typing data/neals.thesis.data"
       subset.data.path = "~/Dropbox (Partners HealthCare)/PAID-UP, PNOIT2 Reactive vs Non-reactive paper/TCR Sequencing/Raw data BL T cell subsets"},
       ns580 = {fig.file.path <- "~/thesis/figures"
       raw.file.path <- "~/data/thesis"},
       stop("I don't recognize your username, type Sys.info() to find out what it is.")
)



# Load in the data --> WGS: i need to back track to understand the rationale of each of these data sets
load(paste(raw.file.path, "parsed.data.baseline.rda", sep = "/"))
load(paste(raw.file.path, "enriched.CDR3s.rda", sep = "/"))
load(paste(raw.file.path, "node.pairs.rda", sep = "/")) # output of network.analysis
load(paste(raw.file.path, "top.nmers.rda", sep = "/"))
load(paste(raw.file.path, "top.disc.fourmers.rda", sep = "/"))
load(paste(raw.file.path, "top.disc.fivemers.rda", sep = "/"))

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
pos.parse <- do.call(rbind, pos.parse)

# Change pos parse to aggregate by CDR3 nucleotide sequence
pos.parse.aggr <- select(pos.parse, c("Read.count", "CDR3.nucleotide.sequence", "CDR3.amino.acid.sequence")) %>%
  aggregate(.~ CDR3.nucleotide.sequence + CDR3.amino.acid.sequence, data = ., sum)

## WGS: ok that length(pos.parse.aggr$CDR3.amino.acid.sequence) > length(unique(pos.parse.aggr$CDR3.amino.acid.sequence))?
## NS: Yes, some unqiue CDR3 AA sequences have multiple nucleotide sequences.  The aggregate function is aggregating read
## counts of rows with the same CDR3 AA AND nucleotide sequence

# Get the CDR3s to look at
net.CDR3s <- unique(c(pairs$from.cdr3, pairs$to.cdr3))

# Set cores for parallel processing
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)
# Need to figure out how many neighbors each node has
neighbors.enriched <- foreach(i = 1:length(net.CDR3s), .combine = c) %dopar%{
  library(stringdist)
  neighbor.fun <- function(x){
  # Determine number of neighbors by homology
  hom <- stringdist(x, pos.parse.aggr$CDR3.amino.acid.sequence, method = "lv")
  # --> WGS: alternative using stringDist::Biostrings (haven't tested yet; note that high score means more homology this will be very sensitive to opening and extension penalty parameters)
  # pairwiseAlignment(pos.parse.aggr$CDR3.amino.acid.sequence, x, substitutionMatrix = "BLOSUM45", gapOpening = 0, gapExtension = 1, scoreOnly = TRUE)
  hom <- hom[hom <= 1]
  # Subtract 1 as to not count itself
  edges <- length(hom) - 1
  
  # Determine number of edges by motifs
  # Determine if CDR3 has a top linear nmer
  if(grepl(paste(top.nmer.df$nmer, sep = "", collapse = "|"), substr(x, 4, nchar(x) - 3))){
    CDR3.vec <- c()
    # Get the nmer
    for(i in 1:length(top.nmer.df$nmer)){
      if (grepl(top.nmer.df$nmer[i], substr(x, 4, nchar(x) - 3))){
        CDR3s <- pos.parse.aggr$CDR3.amino.acid.sequence[grep(top.nmer.df$nmer[i], substr(pos.parse.aggr$CDR3.amino.acid.sequence, 4,
                                                                                          nchar(pos.parse.aggr$CDR3.amino.acid.sequence) - 3))]
        CDR3.vec <- c(CDR3.vec, CDR3s)
      }
    }
    # Get all CDR3s with that nmer
    CDR3.vec <- unique(CDR3.vec)
    
    # Get rid of any with a lev distance of 0 or 1 (already been counted)
    CDR3.vec <- CDR3.vec[stringdist(x, CDR3.vec) > 1]
    
    # Count the edges
    CDR3.count.cont <- nrow(pos.parse.aggr[pos.parse.aggr$CDR3.amino.acid.sequence %in% CDR3.vec,])
    
    # Add to the edges count
    edges <- edges + CDR3.count.cont
    
  }
  # Determine number of edges by discontinuous motifs
  if(nchar(x) >= 10) {
    # Get all discontinuous 4mers
    disc.fourmers <- disc.4mers(data.frame(CDR3 = x, count = 1, stringsAsFactors = FALSE),
                                CDR3.col = "CDR3", count.col = "count")
    # Determine if any are in our list of top discontinuous 4mers
    if(any(disc.fourmers$nmer %in% top.disc.fourmers$nmer)){
      disc.CDR3s <- c()
      # If there is the presence of a disc nmer, loop through the top nmers
      for(i in 1:length(top.disc.fourmers$nmer)){
        # Find which top nmers it has
        if(top.disc.fourmers$nmer[i] %in% disc.fourmers$nmer){
          CDR3.match <- find_disc(top.disc.fourmers$nmer[i], pos.parse.aggr$CDR3.amino.acid.sequence, motif.length = 4)
          
          disc.CDR3s <- c(disc.CDR3s, CDR3.match[!CDR3.match %in% disc.CDR3s])
        }
      }
      # Get rid of CDR3s that have been counted already
      if(exists("CDR3.vec")){
        disc.CDR3s <- disc.CDR3s[!disc.CDR3s %in% CDR3.vec]
        disc.CDR3s <- disc.CDR3s[stringdist(x, disc.CDR3s) > 1]
        # Add the new CDR3s to make sure there aren't duplicates later on when looking at discountinous 5mers
        CDR3.vec <- c(CDR3.vec, disc.CDR3s)
      } else{
        disc.CDR3s <- disc.CDR3s[stringdist(x, disc.CDR3s) > 1]
      }
      # Count them up, subtract one to not count itself
      CDR3.count.disc.4mer <- nrow(pos.parse.aggr[pos.parse.aggr$CDR3.amino.acid.sequence %in% disc.CDR3s,])
      
      # Add to final edges count
      edges <- edges + CDR3.count.disc.4mer
    }
  }
  if(nchar(x) >=11){
    disc.fivemers <- disc.5mers(data.frame(CDR3 = x, count = 1, stringsAsFactors = FALSE),
                                CDR3.col = "CDR3", count.col = "count")
    # Determine if any are in our list of top discontinuous 4mers
    if(any(disc.fivemers$nmer %in% top.disc.fivemers$nmer)){
      disc.CDR3s <- c()
      # If there is the presence of a disc nmer, loop through the top nmers
      for(i in 1:length(top.disc.fivemers$nmer)){
        # Find which top nmers it has
        if(top.disc.fivemers$nmer[i] %in% disc.fivemers$nmer){
          CDR3.match <- find_disc(top.disc.fivemers$nmer[i], pos.parse.aggr$CDR3.amino.acid.sequence, motif.length = 5)
          
          disc.CDR3s <- c(disc.CDR3s, CDR3.match[!CDR3.match %in% disc.CDR3s])
        }
      }
      # Get rid of CDR3s that have been counted already
      if(exists("CDR3.vec")){
        disc.CDR3s <- disc.CDR3s[!disc.CDR3s %in% CDR3.vec]
        disc.CDR3s <- disc.CDR3s[stringdist(x, disc.CDR3s) > 1]
      } else{
        disc.CDR3s <- disc.CDR3s[stringdist(x, disc.CDR3s) > 1]
      }
      # Count them up, subtract one to not count itself
      CDR3.count.disc.5mer <- nrow(pos.parse.aggr[pos.parse.aggr$CDR3.amino.acid.sequence %in% disc.CDR3s,])
      
      # Add to final edges count
      edges <- edges + CDR3.count.disc.5mer
    }
  }
   
  return(edges)
  }
neighbors <- neighbor.fun(net.CDR3s[i]) # ultimately this is to get a count of edges instead of a list glm file, right? so can probalby simplify the derivation of that count by aggregating the information in the glm rather than sort of 're-running' the homology/neighbor type functions
neighbors
}
stopCluster(cl)

neighbor.df <- data.frame(CDR3aa = net.CDR3s,
                          neighbors.enriched = neighbors.enriched)


# Now need to look at number of neighbors in resting CDR3s
neg.parse <- data.parse[grep("neg", names(data.parse))] %>%
  do.call(rbind, .) %>% select(., c("Read.count", "CDR3.nucleotide.sequence", "CDR3.amino.acid.sequence")) %>%
  aggregate(.~ CDR3.nucleotide.sequence + CDR3.amino.acid.sequence, data = ., sum)
  
# Set cores for parallel processing
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

neighbors.neg <- foreach(i = 1:length(net.CDR3s), .combine = c) %dopar%{
  library(stringdist)
  neighbor.fun <- function(x){
    # Determine number of neighbors by homology
    hom <- stringdist(x, neg.parse$CDR3.amino.acid.sequence, method = "lv")
    hom <- hom[hom <= 1]
    # Subtract 1 as to not count itself
    edges <- length(hom)
    
    # Determine number of edges by motifs
    # Determine if CDR3 has a top linear nmer
    if(grepl(paste(top.nmer.df$nmer, sep = "", collapse = "|"), substr(x, 4, nchar(x) - 3))){
      CDR3.vec <- c()
      # Get the nmer
      for(i in 1:length(top.nmer.df$nmer)){
        if (grepl(top.nmer.df$nmer[i], substr(x, 4, nchar(x) - 3))){
          CDR3s <- neg.parse$CDR3.amino.acid.sequence[grep(top.nmer.df$nmer[i], substr(neg.parse$CDR3.amino.acid.sequence, 4,
                                                                                            nchar(neg.parse$CDR3.amino.acid.sequence) - 3))]
          CDR3.vec <- c(CDR3.vec, CDR3s)
        }
      }
      # Get all CDR3s with that nmer
      CDR3.vec <- unique(CDR3.vec)
      
      # Get rid of any with a lev distance of 0 or 1 (alread been counted)
      CDR3.vec <- CDR3.vec[stringdist(x, CDR3.vec) > 1]
      
      # Count the edges
      CDR3.count.cont <- nrow(neg.parse[neg.parse$CDR3.amino.acid.sequence %in% CDR3.vec,])
      
      # Add to the edges count
      edges <- edges + CDR3.count.cont
      
    }

    # Determine number of edges by discontinuous motifs
    if(nchar(x) >= 10) {
      # Get all discontinuous 4mers
      disc.fourmers <- disc.4mers(data.frame(CDR3 = x, count = 1, stringsAsFactors = FALSE),
                                  CDR3.col = "CDR3", count.col = "count")
      # Determine if any are in our list of top discontinuous 4mers
      if(any(disc.fourmers$nmer %in% top.disc.fourmers$nmer)){
        disc.CDR3s <- c()
        # If there is the presence of a disc nmer, loop through the top nmers
        for(i in 1:length(top.disc.fourmers$nmer)){
          # Find which top nmers it has
          if(top.disc.fourmers$nmer[i] %in% disc.fourmers$nmer){
            CDR3.match <- find_disc(top.disc.fourmers$nmer[i], neg.parse$CDR3.amino.acid.sequence, motif.length = 4)
            
            disc.CDR3s <- c(disc.CDR3s, CDR3.match[!CDR3.match %in% disc.CDR3s])
          }
        }
        # Get rid of CDR3s that have been counted already
        if(exists("CDR3.vec")){
          disc.CDR3s <- disc.CDR3s[!disc.CDR3s %in% CDR3.vec]
          disc.CDR3s <- disc.CDR3s[stringdist(x, disc.CDR3s) > 1]
          # Add the new CDR3s to make sure there aren't duplicates later on when looking at discountinous 5mers
          CDR3.vec <- c(CDR3.vec, disc.CDR3s)
        } else{
          disc.CDR3s <- disc.CDR3s[stringdist(x, disc.CDR3s) > 1]
        }
        # Count them up, subtract one to not count itself
        CDR3.count.disc.4mer <- nrow(neg.parse[neg.parse$CDR3.amino.acid.sequence %in% disc.CDR3s,])
        
        # Add to final edges count
        edges <- edges + CDR3.count.disc.4mer
      }
    }
    if(nchar(x) >=11){
      disc.fivemers <- disc.5mers(data.frame(CDR3 = x, count = 1, stringsAsFactors = FALSE),
                                  CDR3.col = "CDR3", count.col = "count")
      # Determine if any are in our list of top discontinuous 4mers
      if(any(disc.fivemers$nmer %in% top.disc.fivemers$nmer)){
        disc.CDR3s <- c()
        # If there is the presence of a disc nmer, loop through the top nmers
        for(i in 1:length(top.disc.fivemers$nmer)){
          # Find which top nmers it has
          if(top.disc.fivemers$nmer[i] %in% disc.fivemers$nmer){
            CDR3.match <- find_disc(top.disc.fivemers$nmer[i], neg.parse$CDR3.amino.acid.sequence, motif.length = 5)
            
            disc.CDR3s <- c(disc.CDR3s, CDR3.match[!CDR3.match %in% disc.CDR3s])
          }
        }
        # Get rid of CDR3s that have been counted already
        if(exists("CDR3.vec")){
          disc.CDR3s <- disc.CDR3s[!disc.CDR3s %in% CDR3.vec]
          disc.CDR3s <- disc.CDR3s[stringdist(x, disc.CDR3s) > 1]
        } else{
          disc.CDR3s <- disc.CDR3s[stringdist(x, disc.CDR3s) > 1]
        }
        # Count them up, subtract one to not count itself
        CDR3.count.disc.5mer <- nrow(neg.parse[neg.parse$CDR3.amino.acid.sequence %in% disc.CDR3s,])
        
        # Add to final edges count
        edges <- edges + CDR3.count.disc.5mer
      }
    }
    return(edges)
  }
  neighbors <- neighbor.fun(net.CDR3s[i]) # searching for neighbors of the CDR3's of interest within the CD154-; could aslo do this with some other control data set(s)
  neighbors
}
stopCluster(cl)
remove(cores);remove(cl)

# Add negative neighbors to dataframe
neighbor.df$neighbors.neg <- neighbors.neg
neighbor.df$neighbors.neg[neighbor.df$neighbors.neg < 0] <- 0

# Perform a G-test looking at neighbors in psCDR3s and CD154-
gtest.df <- data.frame(neighbors.enriched = neighbor.df$neighbors.enriched,
                       neighbors.neg = neighbor.df$neighbors.neg,
                       non.neighbors.enriched = length(net.CDR3s) - neighbor.df$neighbors.enriched,
                       non.nieghbors.neg = length(neg.parse$CDR3.amino.acid.sequence) - neighbor.df$neighbors.neg)
neighbor.df$p.value <- apply(gtest.df, 1, function(x){
  gt <- GTest(matrix(x, ncol = 2))
  return(gt$p.value)
})

# Adjust the P-values with FDR
neighbor.df <- FDR.function(neighbor.df, 0.05)

# Get significant CDR3s and all of their neighbors
core.CDR3s <- as.character(neighbor.df$CDR3aa[neighbor.df$FDRsig == 1]) # more neighbors than expected by chance

## again, for below, i think you could do someting more succint to get these sequences
get.neighbors <- function(x){
  # Determine number of neighbors by homology
  hom <- stringdist(x, pos.parse.aggr$CDR3.amino.acid.sequence, method = "lv")
  neighbors <- pos.parse.aggr$CDR3.amino.acid.sequence[which(hom <= 1)]
  
  # Determine number of edges by motifs
  # Determine if CDR3 has a top linear nmer
  if(grepl(paste(top.nmer.df$nmer, sep = "", collapse = "|"), substr(x, 4, nchar(x) - 3))){
    CDR3.vec <- c()
    # Get the nmer
    for(i in 1:length(top.nmer.df$nmer)){
      if (grepl(top.nmer.df$nmer[i], substr(x, 4, nchar(x) - 3))){
        CDR3s <- pos.parse.aggr$CDR3.amino.acid.sequence[grep(top.nmer.df$nmer[i], substr(pos.parse.aggr$CDR3.amino.acid.sequence, 4,
                                                                                     nchar(pos.parse.aggr$CDR3.amino.acid.sequence) - 3))]
        CDR3.vec <- c(CDR3.vec, CDR3s)
      }
    }
    # Get all CDR3s with that nmer
    CDR3.vec <- unique(CDR3.vec)
    neighbors <- c(neighbors, CDR3.vec)

  }

  # Determine number of edges by discontinuous motifs
  # Get all discontinuous 4mers
  if(nchar(x) >= 10){
    disc.fourmers <- disc.4mers(data.frame(CDR3 = x, count = 1, stringsAsFactors = FALSE),
                                CDR3.col = "CDR3", count.col = "count")
    # Determine if any are in our list of top discontinuous 4mers
    if(any(disc.fourmers$nmer %in% top.disc.fourmers$nmer)){
      disc.CDR3s <- c()
      # If there is the presence of a disc nmer, loop through the top nmers
      for(i in 1:length(top.disc.fourmers$nmer)){
        # Find which top nmers it has
        if(top.disc.fourmers$nmer[i] %in% disc.fourmers$nmer){
          CDR3.match <- find_disc(top.disc.fourmers$nmer[i],pos.parse.aggr$CDR3.amino.acid.sequence, motif.length = 4)
          
          disc.CDR3s <- c(disc.CDR3s, CDR3.match[!CDR3.match %in% disc.CDR3s])
        }
      }
      
      # Combine discontinuous matches to CDR3 vector
      neighbors <- c(neighbors, disc.CDR3s)
    }
  }
  if(nchar(x) >= 11){
    # Get all discontinuous 4mers
    disc.fivemers <- disc.5mers(data.frame(CDR3 = x, count = 1, stringsAsFactors = FALSE),
                                CDR3.col = "CDR3", count.col = "count")
    # Determine if any are in our list of top discontinuous 4mers
    if(any(disc.fivemers$nmer %in% top.disc.fivemers$nmer)){
      disc.CDR3s <- c()
      # If there is the presence of a disc nmer, loop through the top nmers
      for(i in 1:length(top.disc.fivemers$nmer)){
        # Find which top nmers it has
        if(top.disc.fivemers$nmer[i] %in% disc.fivemers$nmer){
          CDR3.match <- find_disc(top.disc.fivemers$nmer[i],pos.parse.aggr$CDR3.amino.acid.sequence, motif.length = 5)
          
          disc.CDR3s <- c(disc.CDR3s, CDR3.match[!CDR3.match %in% disc.CDR3s])
        }
        # Combine discontinuous matches to CDR3 vector
        neighbors <- c(neighbors, disc.CDR3s)
      }
      # Combine discontinuous matches to CDR3 vector
      neighbors <- c(neighbors, disc.CDR3s)
    }
  }
  
  return(unique(neighbors))
}

# Get the CDR3s of neighborhoods of interest
core.neighborhoods <- lapply(core.CDR3s, get.neighbors)

# Set cores for parallel processing
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)
core.neighborhoods <- foreach(i = 1:length(core.CDR3s)) %dopar% {
  library(stringdist)
  data <- get.neighbors(core.CDR3s[i])
  data
}
stopCluster(cl)
remove(cores);remove(cl)
names(core.neighborhoods) <- core.CDR3s

# Save this to be used to probe Teff and Treg compartments
save(core.neighborhoods, file = paste(raw.file.path, "top.neighborhoods.rda", sep = "/"))

