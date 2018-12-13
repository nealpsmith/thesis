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
       raw.file.path <- "C:/Users/nealp/Dropbox (Partners HealthCare)/Projects/PNOIT2-1037/TCRB sequencing and HLA typing data"},
       wayne1 = {fig.file.path <- "~/Dropbox (Personal)/Extension School/Thesis/figures"
       raw.file.path <- "~/Dropbox (Partners HealthCare)/Projects/PNOIT2-1037/TCRB sequencing and HLA typing data"},
       ns580 = {fig.file.path <- "~/thesis/figures"
       raw.file.path <- "~/data/thesis"},
       stop("I don't recognize your username, type Sys.info() to find out what it is.")
)
source(paste(raw.file.path, "neals_tcr_functions.R", sep = "/"))

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


# Function to find CDR3s that are within 0-2 hamming distance of eachother
find_pairs_ham <- function(x, y) {
  res <- dist.to.df(as.dist(stringdistmatrix(x, y,
                                             method = "hamming",
                                             useNames = "strings"))) %>%
    filter(dist == 1 | dist == 2) %>%
    select(-dist)
  colnames(res) <- c("from.cdr3", "to.cdr3")
  res
}

# Find the pairs of CDR3s
# Note that this produces lower number of pairs because it a unique CDR3 from each individual is only counted once
# Earlier, was using same CDR3s with degenerate nucleotide sequences (still could go back to that!)
pairs.ham <- find_pairs_ham(enriched.CDR3s.df$CDR3.amino.acid.sequence, enriched.CDR3s.df$CDR3.amino.acid.sequence)


# Set cores for parallel processing
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

pairs.disc <- foreach(i = 1:length(top.disc.threemers$nmer), .combine = rbind) %dopar% {
  source(paste(raw.file.path, "neals_tcr_functions.R", sep = "/"))
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

pairs.threemers <- foreach(i = 1:length(top.nmers$threemer$nmer), .combine = rbind) %dopar% {
  source(paste(raw.file.path, "neals_tcr_functions.R", sep = "/"))
  data = find_pairs_cont(top.nmers$threemer$nmer[i], enriched.CDR3s.df,
                         CDR3.col = "CDR3.amino.acid.sequence")
  data
}

pairs.fourmers <- foreach(i = 1:length(top.nmers$fourmer$nmer), .combine = rbind) %dopar% {
  source(paste(raw.file.path, "neals_tcr_functions.R", sep = "/"))
  data = find_pairs_cont(top.nmers$fourmer$nmer[i], enriched.CDR3s.df,
                         CDR3.col = "CDR3.amino.acid.sequence")
  data
}

pairs.fivemers <- foreach(i = 1:length(top.nmers$fivemer$nmer), .combine = rbind) %dopar% {
  source(paste(raw.file.path, "neals_tcr_functions.R", sep = "/"))
  data = find_pairs_cont(top.nmers$fivemer$nmer[i], enriched.CDR3s.df,
                         CDR3.col = "CDR3.amino.acid.sequence")
  data
}
stopCluster(cl)


# Add all of the pairs
pairs <- rbind(pairs.ham, pairs.disc, pairs.threemers, pairs.fourmers, pairs.fivemers)

# Try to get edges weighted
layout_graph <- function(graph) {
  set.seed(42)
  # Create a graph from the dataframe
  gg <- graph %>%
    graph_from_data_frame %>%
    simplify
  # Creates clusters: Connects nodes and edges for as long as it can for each CDR3
  cc <- components(gg)
  coords <- gg %>%
    layout_with_mds()
  #layout_with_graphopt(niter = 3000, charge = 0.005)
  data.frame(CDR3aa = names(V(gg)),
             x = coords[,1],
             y = coords[,2],
             stringsAsFactors = F) %>%
    merge(
      data.frame(CDR3aa = names(cc$membership),
                 cid = cc$membership))
}

# Create the graph and also get "frequency" information: number of reads with a CDR3 / total reads
df.mds <- pairs %>%
  layout_graph(.) %>%
  ungroup

# Determine most dominant clusters: Currently clusters with at least 20 CDR3s within hamming of 2 of eachother
dom.clusters <- count(df.mds, vars = cid)$vars[
  count(df.mds, vars = cid)$n >= 10]

ggplot(df.mds, aes(x = x, y = y)) +
  geom_point(data = df.mds[df.mds$cid %in% dom.clusters,], aes(color = as.integer(factor(cid))), alpha = 0.9, shape = 21) +
  geom_point(data = df.mds[!df.mds$cid %in% dom.clusters,], shape = 21, alpha = 0.9) +
  xlab("MDS1") + ylab("MDS2") +
  scale_color_distiller(guide = F, palette = "Set1") +
  scale_size(guide = F) +
  theme_bw() +
  theme(aspect = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
