### TCRNet analysis ###

library(plyr) # load before dplyr to avoid conflicts due to masking
library(dplyr)
library(tcR)
library(ggplot2)
library(GGally)
library(ggseqlogo)
library(Biostrings)
library(msa)
library(grid)
library(gridBase)
library(ggrepel)
library(magrittr)
library(stringi)
library(stringr)
library(stringdist)
source("tcR_functions.R") #includes several modifications to tcR

# Get file path for output figures and raw data files
fig.file.path <- "C:/Users/nealp/Dropbox (Personal)/Extension School/Thesis/figures"
raw.file.path <- "C:/Users/nealp/Dropbox (Partners HealthCare)/Projects/PNOIT2-1037/TCRB sequencing and HLA typing data"

source(paste(raw.file.path, "neals_tcr_functions.R", sep = "/"))

# Load in the data
load(paste(raw.file.path, "neals.thesis.data/parsed.data.baseline.rda", sep = "/"))
load(paste(raw.file.path, "neals.thesis.data/enriched.CDR3s.rda", sep = "/"))

# Adjust format of data to fit TCRNet criteria
# Get just positive data
pos.parse <- data.parse[grep("pos", names(data.parse))]
names(pos.parse) <- substr(names(pos.parse), 0, nchar(names(pos.parse)) - 4)

# # TEST: Make sure names are correct and in order
# names(pos.parse) == names(enriched.CDR3s)

for(i in 1:length(pos.parse)){
  pos.parse[[i]] <- pos.parse[[i]][pos.parse[[i]]$CDR3.amino.acid.sequence %in%
    enriched.CDR3s[[i]]$CDR3.amino.acid.sequence,]
}
pos.parse <- do.call(rbind, pos.parse)

# Get rid of some columns
pos.parse <- pos.parse[,-c(1, 2, 14, 13, 12)]

# Fix the V-gene assaignments
# Going to move V.ties to V.genes for cases with ambiguous assignments
pos.parse$V.gene[pos.parse$V.gene == ""] <- pos.parse$V.ties[pos.parse$V.ties != ""]

pos.parse$V.ties <- NULL

colnames(pos.parse) <- c("count", "frequency", "CDR3nt", "CDR3aa", "V", "J", "D", "Jstart")
pos.parse <- select(pos.parse, c(1:5, 7, 6, 8))

write.csv(pos.parse, "TCRNet.csv")

# Change frequency: Make it the frequency of that CDR3 in the pool of enriched CDR3s
#pos.parse$frequency <- pos.parse$count / sum(pos.parse$count)


# Function to make a distance matrix a dataframe with every row being 2 CDR3s and their distance
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
# Modified from tutorial
find_pairs <- function(x, y) {
  res <- dist.to.df(as.dist(stringdistmatrix(x, y,
                          method = "hamming",
                          useNames = "strings"))) %>%
    filter(dist == 1 | dist == 2 ) %>%
    select(-dist)
  colnames(res) <- c("from.cdr3", "to.cdr3")
  res
}

# Find the pairs of CDR3s
pairs <- find_pairs(pos.parse$CDR3aa, pos.parse$CDR3aa)

#  Create a graphing function
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
# Currently, summing frequencies when CDR3 found across sample: DEF not good idea, will work on soon
df.mds <- pairs %>%
  layout_graph(.) %>%
  ungroup

# Determine most dominant clusters: Currently clusters with at least 20 CDR3s within hamming of 2 of eachother
dom.clusters <- count(df.mds, vars = cid)$vars[
  count(df.mds, vars = cid)$n >=20]

pdf(paste(fig.file.path, "CDR3.clustering.ham.0to2.pdf", sep = "/"), 7, 7)
ggplot(df.mds, aes(x = x, y = y)) +
  geom_point(shape = 21) +
  geom_point(data = df.mds[df.mds$cid %in% dom.clusters,], aes(color = as.integer(factor(cid))), alpha = 0.9) +
  xlab("MDS1") + ylab("MDS2") +
  scale_color_distiller(guide = F, palette = "Set1") +
  scale_size(guide = F) +
  theme_bw() +
  theme(aspect = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
dev.off()

pdf(paste(fig.file.path, "CDR3.clustering.ham.0to2.hist.pdf", sep = "/"), 10, 7)
hist(table(df.mds$cid))
dev.off()

# Okay lets look closer at the dominant clusters
dom.clust.list <- vector("list", length = length(dom.clusters))
k <- 1
for(i in dom.clusters){
  CDR3s <- df.mds$CDR3aa[df.mds$cid == i]
  ham.matrix <- stringdistmatrix(CDR3s, CDR3s, method = "hamming", useNames = "strings")
  ham.matrix <- ifelse(ham.matrix > 2, 0, 1)
  dom.clust.list[[k]] <- ham.matrix
  k <- k+1
}

# Print the graphs out into a pdf
pdf(paste(fig.file.path, "dominant.cluster.networks.pdf", sep = "/"), 7, 7)
for(i in 1:length(dom.clust.list)){
  print(ggnet(dom.clust.list[[i]]) +
          ggtitle(paste("Cluser:", dom.clusters[i], sep = ""))
        )
}
dev.off()

### Look for clustering in random CDR3s of the same sample size ###

# Get just negative data
neg.data <- data.parse[grep("neg", names(data.parse))] %>%
  do.call(rbind, .) %>%
  select(., CDR3.amino.acid.sequence)

# Get just positive data
pos.data <- data.parse[grep("pos", names(data.parse))] %>%
  do.call(rbind, .) %>%
  select(., CDR3.amino.acid.sequence)

# Resample the CD154 negative data in the same amount as enriched data, look at how many edges there are
neg.edges.vec <- vector(mode = "numeric", length = 50) 
for(i in 1:50){
  random.CDR3s <-  sample(neg.data$CDR3.amino.acid.sequence, nrow(pos.parse))
  pairs.control <- find_pairs(random.CDR3s, random.CDR3s)
  neg.edges.vec[i] <- nrow(pairs.control)
}

# Resample unselected CD154 positive data in the same amount as enriched data and look at edges
pos.edges.vec <- vector(mode = "numeric", length = 50)
for(i in 1:50){
  random.CDR3s <-  sample(pos.data$CDR3.amino.acid.sequence, nrow(pos.parse))
  pairs.control <- find_pairs(random.CDR3s, random.CDR3s)
  pos.edges.vec[i] <- nrow(pairs.control)
}


plot.df <- data.frame(variable = c("enriched CDR3s", "Random CD154+", "Random CD154-"), 
                      edges = c(nrow(pairs), median(pos.edges.vec), median(neg.edges.vec)),
                      sd = c(0, sd(pos.edges.vec), sd(neg.edges.vec)))
plot.df$variable <- factor(plot.df$variable, levels = plot.df$variable)

pdf(paste(fig.file.path, "edges.enriched.vs.random.CDR3s.pdf", sep = "/"),7, 7)
ggplot(plot.df, aes(x = variable, y = edges)) + geom_bar(stat = "identity") +
  geom_errorbar(data = plot.df[plot.df$sd > 0,], aes(ymin = edges - sd, 
                                                     ymax = edges + sd), width = 0.2) +
  ggtitle("Edges in enriched vs. randomly sampled CDR3s") +
  xlab("") +
  theme_bw()
dev.off()

# Create graph for control data
df.mds.control <- pairs.control %>%
  layout_graph(.) %>%
  ungroup

# Determine most dominant clusters in a resample of all CD154+
dom.clusters.control <- count(df.mds.control, vars = cid)$vars[
  count(df.mds.control, vars = cid)$n >=20]

pdf(paste(fig.file.path, "clustering.CD154pos.all.pdf", sep = "/"), 7, 7)
ggplot(df.mds.control, aes(x = x, y = y)) +
  geom_point(shape = 21) +
  geom_point(data = df.mds.control[df.mds.control$cid %in% dom.clusters.control,], aes(color = as.integer(factor(cid))), alpha = 0.9) +
  xlab("MDS1") + ylab("MDS2") +
  scale_color_distiller(guide = F, palette = "Set1") +
  scale_size(guide = F) +
  theme_bw() +
  theme(aspect = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
dev.off()

# Make a cluster graph example for random CD154-
random.CDR3s <-  sample(neg.data$CDR3.amino.acid.sequence, nrow(pos.parse))
pairs.control <- find_pairs(random.CDR3s, random.CDR3s)

# Create graph for control data
df.mds.control <- pairs.control %>%
  layout_graph(.) %>%
  ungroup

# Determine most dominant clusters in a resample of all CD154+
dom.clusters.control <- count(df.mds.control, vars = cid)$vars[
  count(df.mds.control, vars = cid)$n >=20]

pdf(paste(fig.file.path, "clustering.CD154neg.pdf", sep = "/"), 7, 7)
ggplot(df.mds.control, aes(x = x, y = y)) +
  geom_point(shape = 21) +
  geom_point(data = df.mds.control[df.mds.control$cid %in% dom.clusters.control,], aes(color = as.integer(factor(cid))), alpha = 0.9) +
  xlab("MDS1") + ylab("MDS2") +
  scale_color_distiller(guide = F, palette = "Set1") +
  scale_size(guide = F) +
  theme_bw() +
  theme(aspect = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
dev.off()

# Want to get sequences of dominant clusters
# Need to "trim" down to dominant clusters
# Trying indise-out approach: Find CDR3 with most edges and look outward to other CDR3s
test <- df.mds[df.mds$cid == "5",]
#test <- dist.to.df(as.dist(stringdistmatrix(test$CDR3aa, test$CDR3aa, method = "hamming", useNames = "strings")))
#test$dist <- ifelse(test$dist > 2, 0, 1)
test <- stringdistmatrix(test$CDR3aa, test$CDR3aa, method = "hamming", useNames = "strings")
test <- ifelse(test > 2, 0, 1)
# Determine number of edges for each CDR3
x <- apply(test, 1, function(x){
  edges <- sum(x) - 1
  return(edges)
})
CDR3.vec <-  names(test[names(x[x == max(x)]),][test[names(x[x==max(x)]),] > 0])
test.2 <- stringdistmatrix(CDR3.vec, CDR3.vec, method = "hamming", useNames = "strings")
test.2 <- ifelse(test.2 > 2, 0, 1)

ggnet2(test.2, label = TRUE) +
  ggtitle(paste("Cluser:", dom.clusters[i], sep = ""))

get_seqs_cid <- function(clust.id) { 
  df.mds %>% 
    filter(cid == clust.id) %>% 
    .$CDR3aa
}

x <- get_seqs_cid("13")



align_seqs <- function(seqs, cons = F) { 
  x <- seqs %>%
    AAStringSet %>% 
    msa(method = "ClustalW")
    if (cons) { 
      return(msaConsensusSequence(.x)) 
    } else { 
        return(x %>% as.matrix %>%
                 melt %>% mutate(seq_id = Var1, base_id = Var2, aa = value) %>%
                 select(-Var1, -Var2, -value) %>% group_by(seq_id) %>% 
                 mutate(seq = paste0(aa[base_id], collapse = "")) %>% 
                 ungroup) 
      }
}

plot_seqgrid <- function(seqs) { 
  seqs %>% align_seqs %>% 
    ggplot(aes(x=base_id, y=seq_id)) + 
    geom_text(aes(label=aa), size = 3) + 
    scale_x_continuous("", breaks = c(), expand = c(0.105, 0)) + 
    theme_logo() + 
    theme(legend.position = 'none')
  }

# plots sequence logo from multiple alignment 
plot_seqlogo <- function(seqs) { 
  seqs %>% align_seqs %>% .$seq %>% 
    unique %>% ggseqlogo + theme(legend.position = 'none')
  }

