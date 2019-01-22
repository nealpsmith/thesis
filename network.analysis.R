library(plyr) # load before dplyr to avoid conflicts due to masking
library(dplyr)
library(ggplot2)
library(GGally)
library(ggseqlogo)
library(Biostrings)
library(ggrepel)
library(magrittr)
library(network)
library(stringdist)
library(igraph)
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

source("neals_tcr_functions.R")
# Load in the data
load(paste(raw.file.path, "parsed.data.baseline.rda", sep = "/")) # all 27 confirmed PN allergic at baseline
load(paste(raw.file.path, "enriched.CDR3s.rda", sep = "/"))
load(paste(raw.file.path, "top.disc.fourmers.rda", sep = "/")) # discontinuous 4mers
load(paste(raw.file.path, "top.disc.fivemers.rda", sep = "/")) # discontinuous 5mers
load(paste(raw.file.path, "top.nmers.rda", sep = "/")) # continuous 3, 4, and 5 mers (list)

# Make a single dataframe for all of the enriched CDR3s
enriched.CDR3s.df <- do.call(rbind, enriched.CDR3s)

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

pairs.disc.4mer <- foreach(i = 1:length(top.disc.fourmers$nmer), .combine = rbind) %dopar% {
  source("neals_tcr_functions.R")
  data = find_pairs_disc(top.disc.fourmers$nmer[i], enriched.CDR3s.df,
                         CDR3.col = "CDR3.amino.acid.sequence", motif.size = 4)
  data
}
stopCluster(cl)
remove(cores); remove(cl)

# Set cores for parallel processing
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

pairs.disc.5mer <- foreach(i = 1:length(top.disc.fivemers$nmer), .combine = rbind) %dopar% {
  source("neals_tcr_functions.R")
  data = find_pairs_disc(top.disc.fivemers$nmer[i], enriched.CDR3s.df,
                         CDR3.col = "CDR3.amino.acid.sequence", motif.size = 5)
  data
}
stopCluster(cl)
remove(cores); remove(cl)


# Set cores for parallel processing
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

# Get pairs for dominant 3mers
pairs.threemers <- foreach(i = 1:length(top.nmers$threemer$nmer), .combine = rbind) %dopar% {
  source("neals_tcr_functions.R") 
  data = find_pairs_cont(top.nmers$threemer$nmer[i], enriched.CDR3s.df,
                         CDR3.col = "CDR3.amino.acid.sequence")
  data
}

# Get pairs for dominant 4mers
pairs.fourmers <- foreach(i = 1:length(top.nmers$fourmer$nmer), .combine = rbind) %dopar% {
  source("neals_tcr_functions.R") 
  data = find_pairs_cont(top.nmers$fourmer$nmer[i], enriched.CDR3s.df,
                         CDR3.col = "CDR3.amino.acid.sequence")
  data
}

# Get pairs for dominant 5mers
pairs.fivemers <- foreach(i = 1:length(top.nmers$fivemer$nmer), .combine = rbind) %dopar% {
  source("neals_tcr_functions.R") 
  data = find_pairs_cont(top.nmers$fivemer$nmer[i], enriched.CDR3s.df,
                         CDR3.col = "CDR3.amino.acid.sequence")
  data
}
stopCluster(cl)
remove(cores); remove(cl)

# Add all of the pairs
pairs <- rbind(pairs.lev, pairs.disc.4mer, pairs.disc.5mer, pairs.threemers, pairs.fourmers, pairs.fivemers)

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
# Need to get correct nucleotide sequence info from the original parsed data set loaded above
pos.parse <- data.parse[grep("pos", names(data.parse))]
names(pos.parse) <- substr(names(pos.parse), 0, nchar(names(pos.parse)) - 4)

# # TEST: Make sure names are correct and in order
# names(pos.parse) == names(enriched.CDR3s) ------> could alternatively encode name matching of some sort

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

# Save the dataframe of pairs of nodes
save(pairs, file = paste(raw.file.path, "node.pairs.rda", sep = "/"))

# Write the pairs to a gml file to examine in cytoscape
pairs.graph <- graph_from_data_frame(pairs)
write_graph(pairs.graph, paste(raw.file.path, "lev1.motif.graph.gml", sep = "/"), format = "gml")


### As a control, Look for the number of edges we see from repeated resampling of random CDR3s ###
# Get just negative data
neg.data <- data.parse[grep("neg", names(data.parse))] %>%
  do.call(rbind, .) %>%
  select(., CDR3.nucleotide.sequence, CDR3.amino.acid.sequence, Read.count) %>%
  apply(., 1, function(y){
  y <- do.call(rbind, replicate(y[["Read.count"]], y, simplify = FALSE))
}) %>% do.call(rbind, .) %>% as.data.frame(.)

# # TEST: Make sure this is the correct sized dataframe
# nrow(neg.data) == sum(unlist(lapply(data.parse[grep("neg", names(data.parse))], function(x) sum(x$Read.count))))

# Get just positive data
pos.data <- data.parse[grep("pos", names(data.parse))] %>%
  do.call(rbind, .) %>%
  select(., CDR3.nucleotide.sequence, CDR3.amino.acid.sequence, Read.count) %>%
  apply(., 1, function(y){
    y <- do.call(rbind, replicate(y[["Read.count"]], y, simplify = FALSE))
  }) %>% do.call(rbind, .) %>% as.data.frame(.)

# # TEST: Make sure this is the correct sized dataframe
# nrow(pos.data) == sum(unlist(lapply(data.parse[grep("pos", names(data.parse))], function(x) sum(x$Read.count))))

# Look in random CD154- CDR3s
neg.edges.vec <- vector(mode = "numeric", length = 50) 
neg.lrg.clust <- vector(mode = "numeric", length = 50)

for(i in 1:50){
  # Get random CDR3s
  indx <- sample.int(nrow(neg.data) ,size = length(enriched.CDR3s.df$CDR3.amino.acid.sequence))
  random.CDR3s <-  neg.data[indx,]
  
  # Change from factors to character vectors (should change this earlier on)
  random.CDR3s$CDR3.nucleotide.sequence <- as.character(random.CDR3s$CDR3.nucleotide.sequence)
  random.CDR3s$CDR3.amino.acid.sequence <- as.character(random.CDR3s$CDR3.amino.acid.sequence)
  # Determine the number of pairs from homology
  pairs.lev.ctrl <- find_pairs_hom(random.CDR3s$CDR3.amino.acid.sequence, random.CDR3s$CDR3.amino.acid.sequence)
  
  ### Determine the number of pairs from nmers ###
  # Set cores for parallel processing
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1)
  registerDoParallel(cl)
  
  # Look for enriched discontinuous 4mers in CD154- CDR3s
  pairs.disc.4mer.ctrl <- foreach(j = 1:length(top.disc.fourmers$nmer), .combine = rbind) %dopar% {
    source("neals_tcr_functions.R")
    data = find_pairs_disc(top.disc.fourmers$nmer[j], random.CDR3s,
                           CDR3.col = "CDR3.amino.acid.sequence", motif.size = 4)
    data
  }
  stopCluster(cl)
  remove(cores); remove(cl)
  
  # Set cores for parallel processing
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1)
  registerDoParallel(cl)
  
  # Look for enriched discontinuous 5mers in CD154- CDR3s
  pairs.disc.5mer.ctrl <- foreach(j = 1:length(top.disc.fivemers$nmer), .combine = rbind) %dopar% {
    source("neals_tcr_functions.R")
    data = find_pairs_disc(top.disc.fivemers$nmer[j], random.CDR3s,
                           CDR3.col = "CDR3.amino.acid.sequence", motif.size = 5)
    data
  }
  stopCluster(cl)
  remove(cores); remove(cl)
  
  # Set cores for parallel processing
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1)
  registerDoParallel(cl)
  
  # Look for enriched 3mers in CD154- CDR3s
  pairs.threemers.ctrl <- foreach(j = 1:length(top.nmers$threemer$nmer), .combine = rbind) %dopar% {
    source("neals_tcr_functions.R") 
    data = find_pairs_cont(top.nmers$threemer$nmer[j], random.CDR3s,
                           CDR3.col = "CDR3.amino.acid.sequence")
    data
  }
  stopCluster(cl)
  remove(cores); remove(cl)
  
  # Set cores for parallel processing
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1)
  registerDoParallel(cl)
  
  # Look for enriched 4mers in CD154- CDR3s
  pairs.fourmers.ctrl <- foreach(j = 1:length(top.nmers$fourmer$nmer), .combine = rbind) %dopar% {
    source("neals_tcr_functions.R") 
    data = find_pairs_cont(top.nmers$fourmer$nmer[j], random.CDR3s,
                           CDR3.col = "CDR3.amino.acid.sequence")
    data
  }
  stopCluster(cl)
  remove(cores); remove(cl)
  
  # Set cores for parallel processing
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1)
  registerDoParallel(cl)
  
  # Look for enriched 5mers in CD154- CDR3s
  pairs.fivemers.ctrl <- foreach(j = 1:length(top.nmers$fivemer$nmer), .combine = rbind) %dopar% {
    source("neals_tcr_functions.R") 
    data = find_pairs_cont(top.nmers$fivemer$nmer[j], random.CDR3s,
                           CDR3.col = "CDR3.amino.acid.sequence")
    data
  }
  stopCluster(cl)
  remove(cores); remove(cl)
  
  # Add all of the pairs
  pairs.ctrl <- rbind(pairs.lev.ctrl, pairs.disc.4mer.ctrl, pairs.disc.5mer.ctrl, pairs.threemers.ctrl,
                      pairs.fourmers.ctrl, pairs.fivemers.ctrl)
  
  # Limit to unique edges
  pairs.ctrl <- unique(pairs.ctrl)
  # Make the columns character vectors
  pairs.ctrl$from.cdr3 <- as.character(pairs.ctrl$from.cdr3)
  pairs.ctrl$to.cdr3 <- as.character(pairs.ctrl$to.cdr3)
  
  # Get rid of any self-edges
  pairs.ctrl <- pairs.ctrl[!pairs.ctrl$from.cdr3 == pairs.ctrl$to.cdr3,]
  
  # Create self-edges for CDR3s with multiple nucleotide sequences in enriched
  # Get all CDR3s in the network
  net.CDR3s.ctrl <- unique(c(pairs.ctrl$from.cdr3, pairs.ctrl$to.cdr3))
  
  # Create a new self edge for every unique nucleotide sequence for the random CDR3s
  self.edges.ctrl <- lapply(net.CDR3s.ctrl, function(x){
    info <- random.CDR3s[random.CDR3s$CDR3.amino.acid.sequence == x,]
    nucleotide.count <- length(unique(info$CDR3.nucleotide.sequence))
    df <- data.frame("from.cdr3" = rep(x, nucleotide.count - 1),
                     "to.cdr3" = rep(x, nucleotide.count - 1))
    
    return(df)
  }) %>% do.call(rbind, .)
  
  pairs.ctrl <- rbind(pairs.ctrl, self.edges.ctrl)
  
  # Add the number of edges from the control CDR3s to a vector
  neg.edges.vec[i] <- nrow(pairs.ctrl)
  
  # Look at "large" clusters
  graph.ctrl <- graph_from_data_frame(pairs.ctrl)
  neg.lrg.clust[i] <- length(clusters(graph.ctrl)$csize[clusters(graph.ctrl)$csize >= 5])
}

# Look in random CD154+ CDR3s
pos.edges.vec <- vector(mode = "numeric", length = 50) 
pos.lrg.clust <- vector(mode = "numeric", length = 50)
for(i in 1:50){
  # Get random CDR3s
  indx <- sample.int(nrow(pos.data) ,size = length(enriched.CDR3s.df$CDR3.amino.acid.sequence))
  random.CDR3s <-  pos.data[indx,]
  
  # Change from factors to character vectors (should change this earlier on)
  random.CDR3s$CDR3.nucleotide.sequence <- as.character(random.CDR3s$CDR3.nucleotide.sequence)
  random.CDR3s$CDR3.amino.acid.sequence <- as.character(random.CDR3s$CDR3.amino.acid.sequence)
  # Determine the number of pairs from homology
  pairs.lev.ctrl <- find_pairs_hom(random.CDR3s$CDR3.amino.acid.sequence, random.CDR3s$CDR3.amino.acid.sequence)
  
  ### Determine the number of pairs from nmers ###
  # Set cores for parallel processing
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1)
  registerDoParallel(cl)
  
  # Look for enriched discontinuous 4mers in CD154- CDR3s
  pairs.disc.4mer.ctrl <- foreach(j = 1:length(top.disc.fourmers$nmer), .combine = rbind) %dopar% {
    source("neals_tcr_functions.R")
    data = find_pairs_disc(top.disc.fourmers$nmer[j], random.CDR3s,
                           CDR3.col = "CDR3.amino.acid.sequence", motif.size = 4)
    data
  }
  stopCluster(cl)
  remove(cores); remove(cl)
  
  # Set cores for parallel processing
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1)
  registerDoParallel(cl)
  
  # Look for enriched discontinuous 5mers in CD154- CDR3s
  pairs.disc.5mer.ctrl <- foreach(j = 1:length(top.disc.fivemers$nmer), .combine = rbind) %dopar% {
    source("neals_tcr_functions.R")
    data = find_pairs_disc(top.disc.fivemers$nmer[j], random.CDR3s,
                           CDR3.col = "CDR3.amino.acid.sequence", motif.size = 5)
    data
  }
  stopCluster(cl)
  remove(cores); remove(cl)
  
  # Set cores for parallel processing
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1)
  registerDoParallel(cl)
  
  # Look for enriched 3mers in CD154- CDR3s
  pairs.threemers.ctrl <- foreach(j = 1:length(top.nmers$threemer$nmer), .combine = rbind) %dopar% {
    source("neals_tcr_functions.R") 
    data = find_pairs_cont(top.nmers$threemer$nmer[j], random.CDR3s,
                           CDR3.col = "CDR3.amino.acid.sequence")
    data
  }
  stopCluster(cl)
  remove(cores); remove(cl)
  
  # Set cores for parallel processing
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1)
  registerDoParallel(cl)
  
  # Look for enriched 4mers in CD154- CDR3s
  pairs.fourmers.ctrl <- foreach(j = 1:length(top.nmers$fourmer$nmer), .combine = rbind) %dopar% {
    source("neals_tcr_functions.R") 
    data = find_pairs_cont(top.nmers$fourmer$nmer[j], random.CDR3s,
                           CDR3.col = "CDR3.amino.acid.sequence")
    data
  }
  stopCluster(cl)
  remove(cores); remove(cl)
  
  # Set cores for parallel processing
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1)
  registerDoParallel(cl)
  
  # Look for enriched 5mers in CD154- CDR3s
  pairs.fivemers.ctrl <- foreach(j = 1:length(top.nmers$fivemer$nmer), .combine = rbind) %dopar% {
    source("neals_tcr_functions.R") 
    data = find_pairs_cont(top.nmers$fivemer$nmer[j], random.CDR3s,
                           CDR3.col = "CDR3.amino.acid.sequence")
    data
  }
  stopCluster(cl)
  remove(cores); remove(cl)
  
  # Add all of the pairs
  pairs.ctrl <- rbind(pairs.lev.ctrl, pairs.disc.4mer.ctrl, pairs.disc.5mer.ctrl, pairs.threemers.ctrl,
                      pairs.fourmers.ctrl, pairs.fivemers.ctrl)
  
  # Limit to unique edges
  pairs.ctrl <- unique(pairs.ctrl)
  # Make the columns character vectors
  pairs.ctrl$from.cdr3 <- as.character(pairs.ctrl$from.cdr3)
  pairs.ctrl$to.cdr3 <- as.character(pairs.ctrl$to.cdr3)
  
  # Get rid of any self-edges
  pairs.ctrl <- pairs.ctrl[!pairs.ctrl$from.cdr3 == pairs.ctrl$to.cdr3,]
  
  # Create self-edges for CDR3s with multiple nucleotide sequences in enriched
  # Get all CDR3s in the network
  net.CDR3s.ctrl <- unique(c(pairs.ctrl$from.cdr3, pairs.ctrl$to.cdr3))
  
  # Create a new self edge for every unique nucleotide sequence for the random CDR3s
  self.edges.ctrl <- lapply(net.CDR3s.ctrl, function(x){
    info <- random.CDR3s[random.CDR3s$CDR3.amino.acid.sequence == x,]
    nucleotide.count <- length(unique(info$CDR3.nucleotide.sequence))
    df <- data.frame("from.cdr3" = rep(x, nucleotide.count - 1),
                     "to.cdr3" = rep(x, nucleotide.count - 1))
    
    return(df)
  }) %>% do.call(rbind, .)
  
  pairs.ctrl <- rbind(pairs.ctrl, self.edges.ctrl)
  
  # Add the number of edges from the control CDR3s to a vector
  pos.edges.vec[i] <- nrow(pairs.ctrl)
  
  # Look at "large" clusters
  graph.ctrl <- graph_from_data_frame(pairs.ctrl)
  pos.lrg.clust[i] <- length(clusters(graph.ctrl)$csize[clusters(graph.ctrl)$csize >= 5])
  
}

# Create a dataframe to plot 
plot.df <- data.frame(variable = c("Random CD154-", "Random CD154+", "psCDR3s"), 
                      edges = c(median(neg.edges.vec), median(pos.edges.vec), nrow(pairs)),
                      sd = c(sd(neg.edges.vec), sd(pos.edges.vec), 0))
plot.df$variable <- factor(plot.df$variable, levels = plot.df$variable)

# Plot the median number of edges created by homology
pdf(paste(fig.file.path, "edges.psCDR3s.vs.random.CDR3s.pdf", sep = "/"),7, 7)
ggplot(plot.df, aes(x = variable, y = edges)) + geom_bar(stat = "identity") +
  geom_errorbar(data = plot.df[plot.df$sd > 0,], aes(ymin = edges - sd, 
                                                     ymax = edges + sd), width = 0.2) +
  ggtitle("Edges in psCDR3s vs. randomly sampled CDR3s") +
  xlab("") + 
  theme_bw()
dev.off()


x <- data.frame(size = x$csize)
ggplot(x, aes(x = size)) + geom_histogram(bins = 100)


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


