library(plyr) # load before dplyr to avoid conflicts due to masking
library(dplyr)
library(tcR)
library(ggplot2)
library(grid)
library(gridBase)
library(ggrepel)
library(magrittr)
library(data.table)
library(igraph)
library(gplots)
source("tcR_functions.R") #includes several modifications to tcR
source("neals_tcr_functions.R")

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



# Load in the data
load(paste(raw.file.path, "parsed.data.baseline.rda", sep = "/"))
load(paste(raw.file.path, "enriched.CDR3s.rda", sep = "/"))
load(paste(raw.file.path, "top.neighborhoods.rda", sep = "/"))
load(paste(raw.file.path, "node.pairs.rda", sep = "/"))

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

# Need to get some basic stats about the Teff/Treg data
subset.stats <- data.frame("tot.teff.CDR3s" = unlist(lapply(teff.parse, function(x) length(unique(x$CDR3.amino.acid.sequence)))),
                "tot.treg.CDR3s" = unlist(lapply(treg.parse, function(x) length(unique(x$CDR3.amino.acid.sequence)))),
                "tot.teff.reads" = unlist(lapply(teff.parse, function(x) sum(x$Read.count))),
                "tot.treg.reads" = unlist(lapply(treg.parse, function(x) sum(x$Read.count))),
                "psCDR3.teff" = unlist(lapply(teff.psCDR3, function(x) length(unique(x$CDR3.amino.acid.sequence)))),
                "psCDR3.treg" = unlist(lapply(treg.psCDR3, function(x) length(unique(x$CDR3.amino.acid.sequence)))),
                "psCDR3.teff.count" = unlist(lapply(teff.psCDR3, function(x) sum(x$Read.count))),
                "psCDR3.treg.count" = unlist(lapply(treg.psCDR3, function(x) sum(x$Read.count))))

subset.stats$psCDR3.perc.teff <- (subset.stats$psCDR3.teff / subset.stats$tot.teff.CDR3s) * 100
subset.stats$psCDR3.perc.treg <- (subset.stats$psCDR3.treg / subset.stats$tot.treg.CDR3s) * 100
subset.stats$psCDR3.perc.teff.count <- (subset.stats$psCDR3.teff.count / subset.stats$tot.teff.reads) * 100
subset.stats$psCDR3.perc.treg.count <- (subset.stats$psCDR3.treg.count / subset.stats$tot.treg.reads) * 100

plot.df <- melt(subset.stats[c("psCDR3.perc.teff", "psCDR3.perc.treg", "psCDR3.perc.teff.count", "psCDR3.perc.treg.count")])

plot.df$facet <- ifelse(grepl("count", plot.df$variable), "read.count", "unique.CDR3s")

pdf(paste(fig.file.path, "psCDR3.treg.teff.percents.pdf", sep = "/"), 10, 7)
ggplot(plot.df, aes(x = variable, y = value)) + geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1)) +
  ggtitle("psCDR3 Teff/Treg percents") +
  scale_x_discrete(labels = c("Teff", "Treg")) + 
  xlab("") + ylab("Percent") +
  facet_wrap(~facet, scales = "free", labeller = as_labeller(c(read.count = "Read count", unique.CDR3s = "Unique CDR3s"))) +
  theme_bw(base_size = 15)
dev.off()

plot.df <- melt(subset.stats[c("tot.teff.CDR3s", "tot.treg.CDR3s", "tot.teff.reads", "tot.treg.reads")])
plot.df$facet <- ifelse(grepl("reads", plot.df$variable), "read.count", "unique.CDR3s")

pdf(paste(fig.file.path, "teff.treg.reads.pdf", sep = "/"), 10, 7)
ggplot(plot.df, aes(x = variable, y = value)) + geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1)) +
  ggtitle("subset sequences") +
  scale_x_discrete(labels = c("Teff", "Treg")) + 
  xlab("") + ylab("sequences") +
  facet_wrap(~facet, scales = "free", labeller = as_labeller(c(read.count = "Read count", unique.CDR3s = "Unique CDR3s"))) +
  theme_bw(base_size = 15)
dev.off()

write.csv(subset.stats, paste(fig.file.path, "subset.stats.csv", sep = "/"))


### Look at psCDR3s in Teff and Treg data ###

teff.psCDR3s.prop <- lapply(teff.parse, function(x){
  psCDR3.prop <- lapply(unique(enriched.CDR3s.df$CDR3.amino.acid.sequence), function(y){
    psCDR3 <- x[x$CDR3.amino.acid.sequence %in% as.character(y),]
    prop <- sum(psCDR3$Read.count) / sum(x$Read.count)
    return(prop)
  })
  return(psCDR3.prop)
})

# Create a dataframe with this data
teff.psCDR3s.prop <- as.data.frame(matrix(unlist(teff.psCDR3s.prop), nrow = length(unlist(teff.psCDR3s.prop[1]))))
colnames(teff.psCDR3s.prop) <- names(teff.parse)
teff.psCDR3s.prop$psCDR3 <- unique(enriched.CDR3s.df$CDR3.amino.acid.sequence)
teff.psCDR3s.prop <- select(teff.psCDR3s.prop, psCDR3, everything())

# Determine the number of Teff samples that have a CDR3
teff.psCDR3s.prop$samp.count.teff <- apply(teff.psCDR3s.prop[,2:ncol(teff.psCDR3s.prop)], 1, function(x) length(x[x != 0]))


treg.psCDR3s.prop <- lapply(treg.parse, function(x){
  psCDR3.prop <- lapply(unique(enriched.CDR3s.df$CDR3.amino.acid.sequence), function(y){
    psCDR3 <- x[x$CDR3.amino.acid.sequence %in% as.character(y),]
    prop <- sum(psCDR3$Read.count) / sum(x$Read.count)
    return(prop)
  })
  return(psCDR3.prop)
})

# Create a dataframe with this data
treg.psCDR3s.prop <- as.data.frame(matrix(unlist(treg.psCDR3s.prop), nrow = length(unlist(treg.psCDR3s.prop[1]))))
colnames(treg.psCDR3s.prop) <- names(treg.parse)
treg.psCDR3s.prop$psCDR3 <- unique(enriched.CDR3s.df$CDR3.amino.acid.sequence)
treg.psCDR3s.prop <- select(treg.psCDR3s.prop, psCDR3, everything())

# Determine the number of Treg samples that have a CDR3 
treg.psCDR3s.prop$samp.count.treg <- apply(treg.psCDR3s.prop[,2:ncol(treg.psCDR3s.prop)], 1, function(x) length(x[x != 0]))


# Join the Teff and Treg psCDR3 data
psCDR3.prop.df <- left_join(teff.psCDR3s.prop, treg.psCDR3s.prop, by = "psCDR3") %>%
  select(c("psCDR3", colnames(.)[grep("aTeff", colnames(.))],
           colnames(.)[grep("iTreg", colnames(.))],
           "samp.count.teff", "samp.count.treg"))

# Get only public psCDR3s among single compartment
psCDR3.prop.df.pub <- psCDR3.prop.df[psCDR3.prop.df$samp.count.teff > 1 | psCDR3.prop.df$samp.count.treg > 1,] %>%
  arrange(samp.count.teff > 0 & samp.count.treg == 0) %>% arrange(samp.count.treg > 0 & samp.count.teff > 0) %>%
  arrange(samp.count.teff == 0 & samp.count.treg > 0)



# Get the rownames of the heatmap
rownames <- psCDR3.prop.df.pub[,1]
# Create a matrix for the heatmap
mat_data <- data.matrix(psCDR3.prop.df.pub[,2:(ncol(psCDR3.prop.df) - 2)])  
rownames(mat_data) <- rownames 
# Create a color palette for the plot
my_palette <- colorRampPalette(c("black", "black", "red"))(n = 1000)

pdf(paste(fig.file.path, "public.psCDR3.heatmap.in.subsets.pdf", sep = "/"), 10, 5)
heatmap.2(mat_data, scale = "row", tracecol=NA, Rowv = F, Colv = F, dendrogram = "none",
          col = my_palette)
dev.off()

### Compare the Teff and Treg psCDR3s by BLAST algorithm ###
# Combine the Teff/Treg psCDR3s
teff.psCDR3.seqs <- lapply(teff.psCDR3, function(x){
  select(x, CDR3.amino.acid.sequence)
}) %>% do.call(rbind, .) %>% .$CDR3.amino.acid.sequence %>% unique()

treg.psCDR3.seqs <- lapply(treg.psCDR3, function(x){
  select(x, CDR3.amino.acid.sequence)
}) %>% do.call(rbind, .) %>% .$CDR3.amino.acid.sequence %>% unique()

subset.psCDR3.seqs <- c(teff.psCDR3.seqs, treg.psCDR3.seqs) %>% unique()

# Determine their scores by BLAST algorithm using BLOSUM62 matrix
blos.mtx <- matrix(nrow = length(subset.psCDR3.seqs), ncol = length(subset.psCDR3.seqs)) %>%
  `rownames<-`(., subset.psCDR3.seqs)

for(i in 1:length(subset.psCDR3.seqs)){
  blos.mtx[,i] <- pairwiseAlignment(subset.psCDR3.seqs, subset.psCDR3.seqs[i], substitutionMatrix = "BLOSUM62")@score
  
}

# Get the order of the sequences
blos.seq.order <- hclust(dist(blos.mtx), method = "ward.D")$order

blos.clust.df <- data.frame(psCDR3 = subset.psCDR3.seqs[blos.seq.order]) %>%
  left_join(psCDR3.prop.df, by = "psCDR3")
blos.clust.df$mean.teff <- rowMeans(select(blos.clust.df, grep("aTeff", colnames(blos.clust.df))))
blos.clust.df$mean.treg <- rowMeans(select(blos.clust.df, grep("iTreg", colnames(blos.clust.df))))

rownames <- rownames(blos.clust.df$psCDR3)
mat_data <- data.matrix(select(blos.clust.df, grep("mean", colnames(blos.clust.df))))  
rownames(mat_data) <- rownames 
# Create a color palette for the plot
my_palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 1000)

heatmap.2(mat_data, scale = "row", tracecol=NA, Rowv = F, Colv = F, dendrogram = "none",
          col = my_palette)




# Look for CDR3s of top neighborhoods in Teff
teff.neighborhoods <- lapply(teff.parse, function(x){
  # Iterate through all of the neighborhoods
  neighborhood.prop <- lapply(core.neighborhoods, function(y){
    CDR3s <- x[x$CDR3.amino.acid.sequence %in% as.character(y),]
    prop <- sum(CDR3s$Read.count) / sum(x$Read.count)
    return(prop)
  })
  return(neighborhood.prop)
})

# Create a dataframe with this data
teff.neighborhoods <- as.data.frame(matrix(unlist(teff.neighborhoods), nrow = length(unlist(teff.neighborhoods[1]))))
colnames(teff.neighborhoods) <- names(teff.parse)
teff.neighborhoods$neighborhood <- as.character(seq(1:length(core.neighborhoods)))
teff.neighborhoods <- select(teff.neighborhoods, neighborhood, everything())

# Determine the number of Teff samples that have a CDR3 in a given neighborhood
teff.neighborhoods$samp.count.teff <- apply(teff.neighborhoods[,2:ncol(teff.neighborhoods)], 1, function(x) length(x[x != 0]))

# Look for CDR3s of top neighborhoods in Teff/Treg data
treg.neighborhoods <- lapply(treg.parse, function(x){
  # Iterate through all of the neighborhoods
  neighborhood.prop <- lapply(core.neighborhoods, function(y){
    CDR3s <- x[x$CDR3.amino.acid.sequence %in% as.character(y),]
    prop <- sum(CDR3s$Read.count) / sum(x$Read.count)
    return(prop)
  })
  return(neighborhood.prop)
})
# Create a dataframe with this data
treg.neighborhoods <- as.data.frame(matrix(unlist(treg.neighborhoods), nrow = length(unlist(treg.neighborhoods[1]))))
colnames(treg.neighborhoods) <- names(treg.parse)
treg.neighborhoods$neighborhood <- as.character(seq(1:length(core.neighborhoods)))
treg.neighborhoods <- select(treg.neighborhoods, neighborhood, everything())

# Determine the number of Treg samples that have a CDR3 in a given neighborhood
treg.neighborhoods$samp.count.treg <- apply(treg.neighborhoods[,2:ncol(treg.neighborhoods)], 1, function(x) length(x[x != 0]))

# Join the Teff and Treg neighborhood
neighborhood.df <- left_join(teff.neighborhoods, treg.neighborhoods, by = "neighborhood") %>%
  select(c("neighborhood", colnames(.)[grep("aTeff", colnames(.))],
           colnames(.)[grep("iTreg", colnames(.))],
           "samp.count.teff", "samp.count.treg")) %>%
  # First, arrange bt presence in Teff
  arrange(samp.count.teff)

# Get rid of any neighborhoods not present in subset data
neighborhood.df <- neighborhood.df[neighborhood.df$samp.count.teff > 0 | neighborhood.df$samp.count.treg > 0,]

# Arrange the data frame such that those entirely specific to a single compartment are together on the heatmap
neighborhood.df <- arrange(neighborhood.df, samp.count.treg)

# Get the rownames of the heatmap
rownames <- neighborhood.df[,1]
# Create a matrix for the heatmap
mat_data <- data.matrix(neighborhood.df[,2:(ncol(neighborhood.df) - 2)])  
rownames(mat_data) <- rownames 
# Create a color palette for the plot
my_palette <- colorRampPalette(c("black", "black", "red"))(n = 1000)

# Create a heatmap where each row is scaled
pdf(paste(fig.file.path, "neighborhood.heatmap.in.subsets.pdf", sep = "/"), 10, 5)
heatmap.2(mat_data, scale = "row", tracecol=NA, dendrogram = "none", Colv = FALSE, Rowv = FALSE,
          col = my_palette)
dev.off()


# Make a graph with the Teff/Treg information
pairs.graph <- graph_from_data_frame(pairs)
# Add in node attributes based on Teff or Treg presence
pairs.graph <- igraph::as_data_frame(pairs.graph, what = "both")

# Want to make less binary, determine if any are "skewed" to a subset
pairs.graph$vertices$subset.prop <- unlist(lapply(pairs.graph$vertices$name, function(x){
  teff.prop.mean <- mean(unlist(lapply(teff.parse, function(y){
    CDR3s <- y[y$CDR3.amino.acid.sequence == x,]
    prop <- sum(CDR3s$Read.count) / sum(y$Read.count)
    return(prop)
  })))
  treg.prop.mean <- mean(unlist(lapply(treg.parse, function(y){
    CDR3s <- y[y$CDR3.amino.acid.sequence == x,]
    prop <- sum(CDR3s$Read.count) / sum(y$Read.count)
    return(prop)
  })))
  
  # Log adjust to make it a more reasonable scale
  prop.log <- log(teff.prop.mean / treg.prop.mean)
  if(is.nan(prop.log)){
    return(0)
  } else if(prop.log == "Inf") {
    return(log(50))
  } else if (prop.log == "-Inf"){
    return(log(0.02))
  } else{
    return(prop.log)
  }

}))

pairs.graph$vertices$subset.prop <- pairs.graph$vertices$subset.prop + 10.01

# Re-make the graph with the newly added attributes
pairs.graph <- graph_from_data_frame(pairs.graph$edges,
                                     directed = F,
                                     vertices = pairs.graph$vertices)

# Write the graph to be analyzed in cytoscape
write_graph(pairs.graph,paste(raw.file.path, "lev1.motif.graph.with.compartment.data.gml", sep = "/"), format = "gml")


### Examine motifs in Teff/Treg subsets ###
load(paste(raw.file.path, "top.nmers.rda", sep = "/"))
load(paste(raw.file.path, "top.disc.fourmers.rda", sep = "/"))

# Get all continuous nmers into a single dataframe
top.cont.nmer.df <- do.call(rbind, top.nmers)

# Search for CDR3s with the nmers in the subsets
# First look in the Teff
cont.motif.prop.Teff <- lapply(teff.parse, function(x){
  # Search for the nmer
  search <- sapply(top.cont.nmer.df$nmer, function(y){
    # Determine if the nmer is present
    lookup <- nmer.search(y, x, CDR3.col = 6, count.col = 3)
    return(lookup$clone.count)
  })
  names(search) <- top.cont.nmer.df$nmer
  search <- lapply(search, function(z) z / sum(x$Read.count))
  return(search)
}) 

cont.motif.prop.Teff <- data.frame(matrix(unlist(cont.motif.prop.Teff), nrow = length(top.cont.nmer.df$nmer)),
                                   nmer = top.cont.nmer.df$nmer) %>%
  select(., "nmer", 1:(ncol(.) - 1)) %>%
  set_colnames(c("nmer", names(teff.parse)))

# Determine the number of Teff samples that have a CDR3 with the 
cont.motif.prop.Teff$samp.count.teff <- apply(cont.motif.prop.Teff[,2:ncol(cont.motif.prop.Teff)], 1, function(x) length(x[x != 0]))


cont.motif.prop.Treg <- lapply(treg.parse, function(x){
  # Search for the nmer
  search <- sapply(top.cont.nmer.df$nmer, function(y){
    # Determine if the nmer is present
    lookup <- nmer.search(y, x, CDR3.col = 6, count.col = 3)
    return(lookup$clone.count)
  })
  names(search) <- top.cont.nmer.df$nmer
  search <- lapply(search, function(z) z / sum(x$Read.count))
  return(search)
}) 

cont.motif.prop.Treg <- data.frame(matrix(unlist(cont.motif.prop.Treg), nrow = length(top.cont.nmer.df$nmer)),
                                   nmer = top.cont.nmer.df$nmer) %>%
  select(., "nmer", 1:(ncol(.) - 1)) %>%
  set_colnames(c("nmer", names(treg.parse)))

# Determine the number of Treg samples that have a CDR3 in a given neighborhood
cont.motif.prop.Treg$samp.count.treg <- apply(cont.motif.prop.Treg[,2:ncol(cont.motif.prop.Treg)], 1, function(x) length(x[x != 0]))

# Join the Teff and Treg data
motif.prop.df <- left_join(cont.motif.prop.Teff, cont.motif.prop.Treg, by = "nmer") %>%
  arrange(samp.count.teff)
# Arrange the data frame such that those entirely specific to a single compartment are together on the heatmap
 motif.prop.df <- arrange(motif.prop.df, samp.count.treg) %>%
   select("nmer", grep("aTeff", colnames(.)), grep("iTreg", colnames(.)), "samp.count.teff", "samp.count.treg")


# Get the rownames of the heatmap
rownames <- motif.prop.df[,1]
# Create a matrix for the heatmap
mat_data <- data.matrix(motif.prop.df[,2:(ncol(motif.prop.df) - 2)])  
rownames(mat_data) <- rownames 
# Create a color palette for the plot
my_palette <- colorRampPalette(c("black", "black", "red"))(n = 1000)

# Create a heatmap where each row is scaled
pdf(paste(fig.file.path, "neighborhood.heatmap.in.subsets.pdf", sep = "/"), 10, 5)
heatmap.2(mat_data, scale = "row", tracecol=NA, dendrogram = "none", Rowv = FALSE, Colv = FALSE,
          col = my_palette)
dev.off()


