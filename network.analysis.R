library(plyr) # load before dplyr to avoid conflicts due to masking
library(dplyr)
library(tcR)
library(ggplot2)
library(GGally)
library(ggseqlogo)
library(Biostrings)
#library(msa)
library(ggrepel)
library(magrittr)
library(network)
library(stringdist)
library(doParallel)
library(foreach)

switch(Sys.info()[['user']],
       nealp = {fig.file.path <- "C:/Users/nealp/Dropbox (Personal)/Extension School/Thesis/figures"
       raw.file.path <- "C:/Users/nealp/Dropbox (Partners HealthCare)/Projects/PNOIT2-1037/TCRB sequencing and HLA typing data/neals.thesis.data"},
       wayne1 = {fig.file.path <- "~/Dropbox (Personal)/Extension School/Thesis/figures"
       raw.file.path <- "~/Dropbox (Partners HealthCare)/Projects/PNOIT2-1037/TCRB sequencing and HLA typing data/neals.thesis.data"},
       ns580 = {fig.file.path <- "~/thesis/figures"
       raw.file.path <- "~/data/thesis"},
       stop("I don't recognize your username, type Sys.info() to find out what it is.")
)

source("~/data/neals_tcr_functions.R")

# Load in the data
load(paste(raw.file.path, "parsed.data.baseline.rda", sep = "/"))
load(paste(raw.file.path, "enriched.CDR3s.rda", sep = "/"))
load(paste(raw.file.path, "top.disc.threemers.rda", sep = "/"))
load(paste(raw.file.path, "top.nmers.rda", sep = "/"))

# Make a single dataframe for all of the enriched CDR3s
enriched.CDR3s.df <- do.call(rbind, enriched.CDR3s)

# Needed for next function
dist.to.df <- function(inDist) {
  if (class(inDist) != "dist") stop("wrong input type")
  A <- attr(inDist, "Size")
  B <- if (is.null(attr(inDist, "Labels"))) sequence(A) else attr(inDist, "Labels")
  if (isTRUE(attr(inDist, "Diag"))) attr(inDist, "Diag") <- FALSE
  if (isTRUE(attr(inDist, "Upper"))) attr(inDist, "Upper") <- FALSE
  data.frame(
    row = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
    col = rep(B[-length(B)], (length(B)-1):1),
    dist = as.vector(inDist))
}


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

# Find the pairs of CDR3s
# Note that this produces lower number of pairs because it a unique CDR3 from each individual is only counted once
# Earlier, was using same CDR3s with degenerate nucleotide sequences (still could go back to that!)
pairs.lev <- find_pairs_hom(enriched.CDR3s.df$CDR3.amino.acid.sequence, enriched.CDR3s.df$CDR3.amino.acid.sequence)


# Set cores for parallel processing
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

pairs.disc <- foreach(i = 1:length(top.disc.threemers$nmer), .combine = rbind) %dopar% {
  source("~/data/neals_tcr_functions.R")
  data = find_pairs_disc(top.disc.threemers$nmer[i], enriched.CDR3s.df,
                         CDR3.col = "CDR3.amino.acid.sequence")
  data
}
stopCluster(cl)

# Get pairs based on continuous nmers
find_pairs_cont <- function(nmer, df, CDR3.col){
  
  # Iterate across the CDR3s looking for which ones have an nmer
  CDR3s <- df[[CDR3.col]][grepl(nmer, substr(df[[CDR3.col]], 4, nchar(df[[CDR3.col]]) - 3))]
  pairs <- combn(CDR3s, 2)
  df <- data.frame("from.cdr3" = pairs[1,],
                   "to.cdr3" = pairs[2,])
  return(df)
}

# Set cores for parallel processing
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

# Get pairs for dominant 3mers
pairs.threemers <- foreach(i = 1:length(top.nmers$threemer$nmer), .combine = rbind) %dopar% {
  source("~/data/neals_tcr_functions.R")
  data = find_pairs_cont(top.nmers$threemer$nmer[i], enriched.CDR3s.df,
                         CDR3.col = "CDR3.amino.acid.sequence")
  data
}

# Get pairs for dominant 4mers
pairs.fourmers <- foreach(i = 1:length(top.nmers$fourmer$nmer), .combine = rbind) %dopar% {
  source("~/data/neals_tcr_functions.R")
  data = find_pairs_cont(top.nmers$fourmer$nmer[i], enriched.CDR3s.df,
                         CDR3.col = "CDR3.amino.acid.sequence")
  data
}

# Get pairs for dominant 5mers
pairs.fivemers <- foreach(i = 1:length(top.nmers$fivemer$nmer), .combine = rbind) %dopar% {
  source("~/data/neals_tcr_functions.R")
  data = find_pairs_cont(top.nmers$fivemer$nmer[i], enriched.CDR3s.df,
                         CDR3.col = "CDR3.amino.acid.sequence")
  data
}
stopCluster(cl)


# Add all of the pairs
pairs <- rbind(pairs.lev, pairs.disc, pairs.threemers, pairs.fourmers, pairs.fivemers)

# Limit to unique edges
pairs <- unique(pairs)
# Make the columns character vectors
pairs$from.cdr3 <- as.character(pairs$from.cdr3)
pairs$to.cdr3 <- as.character(pairs$to.cdr3)

# Get rid of any self-edges
pairs <- pairs[!pairs$from.cdr3 == pairs$to.cdr3,]

# Create self-edges for CDR3s with multiple nucleotide sequences in enriched
# Get all CDR3s in the network
net.CDR3s <- unique(c(pairs$from.cdr3, pairs$to.cdr3))

# Look for multiple nucleotide sequences for each CDR3
# Need to get correct nucleotide sequence info
pos.parse <- data.parse[grep("pos", names(data.parse))]
names(pos.parse) <- substr(names(pos.parse), 0, nchar(names(pos.parse)) - 4)

# # TEST: Make sure names are correct and in order
# names(pos.parse) == names(enriched.CDR3s)

# Make a single dataframe that has nucleotide info for everyones enriched CDR3s
for(i in 1:length(pos.parse)){
  pos.parse[[i]] <- pos.parse[[i]][pos.parse[[i]]$CDR3.amino.acid.sequence %in%
                                     enriched.CDR3s[[i]]$CDR3.amino.acid.sequence,]
}
pos.parse <- do.call(rbind, pos.parse)

# Create a new self edge for every unique nucleotide sequence for an enriched CDR3
self.edges <- lapply(net.CDR3s, function(x){
  info <- pos.parse[pos.parse$CDR3.amino.acid.sequence == x,]
  nucleotide.count <- length(unique(info$CDR3.nucleotide.sequence))
  df <- data.frame("from.cdr3" = rep(x, nucleotide.count - 1),
                   "to.cdr3" = rep(x, nucleotide.count - 1))
  
  return(df)
}) %>% do.call(rbind, .)

# Add self-edges to pair list
pairs <- rbind(pairs, self.edges)

# Write the pairs to a gml file to examine in cytoscape
pairs.graph <- graph_from_data_frame(pairs)
write_graph(pairs.graph, "lev1.motif.graph.gml", format = "gml")


### As a control, Look for the number of edges we see from repeated resampling of resting CDR3s ###
# Get just negative data
neg.data <- data.parse[grep("neg", names(data.parse))] %>%
  do.call(rbind, .) %>%
  select(., CDR3.amino.acid.sequence)

# Get just positive data
pos.data <- data.parse[grep("pos", names(data.parse))] %>%
  do.call(rbind, .) %>%
  select(., CDR3.amino.acid.sequence)

neg.edges.vec <- vector(mode = "numeric", length = 50) 
for(i in 1:50){
  random.CDR3s <-  sample(neg.data$CDR3.amino.acid.sequence, length(unique(enriched.CDR3s.df$CDR3.amino.acid.sequence)))
  pairs.control <- find_pairs_hom(random.CDR3s, random.CDR3s)
  neg.edges.vec[i] <- nrow(pairs.control)
}

# Resample unselected CD154 positive data in the same amount as enriched data and look at edges
pos.edges.vec <- vector(mode = "numeric", length = 50)
for(i in 1:50){
  random.CDR3s <-  sample(pos.data$CDR3.amino.acid.sequence, length(unique(enriched.CDR3s.df$CDR3.amino.acid.sequence)))
  pairs.control <- find_pairs_hom(random.CDR3s, random.CDR3s)
  pos.edges.vec[i] <- nrow(pairs.control)
}

# Create a dataframe to plot 
plot.df <- data.frame(variable = c("enriched CDR3s", "Random CD154+", "Random CD154-"), 
                      edges = c(nrow(pairs.lev), median(pos.edges.vec), median(neg.edges.vec)),
                      sd = c(0, sd(pos.edges.vec), sd(neg.edges.vec)))
plot.df$variable <- factor(plot.df$variable, levels = plot.df$variable)

# Plot the median number of edges created by homology
pdf(paste(fig.file.path, "lev.edges.enriched.vs.random.CDR3s.pdf", sep = "/"),7, 7)
ggplot(plot.df, aes(x = variable, y = edges)) + geom_bar(stat = "identity") +
  geom_errorbar(data = plot.df[plot.df$sd > 0,], aes(ymin = edges - sd, 
                                                     ymax = edges + sd), width = 0.2) +
  ggtitle("Edges in enriched vs. randomly sampled CDR3s") +
  xlab("") +
  theme_bw()
dev.off()


### Want to take clusters and look for any HLA association ###
# Load HLA data 
HLA.info <- read.csv(paste(raw.file.path, "HLA.MHCII.alleles.csv", sep = "/"))

# Format IDs so they match other
HLA.info$id[HLA.info$id < 10] <- paste("0", HLA.info$id[HLA.info$id < 10], sep = "")
HLA.info$id[HLA.info$id == "06"] <- "19"

# Get CDR3s from a dominant cluster
x <- enriched.CDR3s.df[grep("", enriched.CDR3s.df$CDR3.amino.acid.sequence),] %>%
  select(c("id", "CDR3.amino.acid.sequence")) %>%
  left_join(HLA.info, by = "id")


