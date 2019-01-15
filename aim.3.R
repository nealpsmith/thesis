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
### Note: Unique CDR3s here is by nucelotide, not AA ###
subset.stats <- data.frame("tot.teff.CDR3s" = unlist(lapply(teff.parse, nrow)),
                "tot.treg.CDR3s" = unlist(lapply(treg.parse, nrow)),
                "tot.teff.reads" = unlist(lapply(teff.parse, function(x) sum(x$Read.count))),
                "tot.treg.reads" = unlist(lapply(treg.parse, function(x) sum(x$Read.count))),
                "psCDR3.teff" = unlist(lapply(teff.psCDR3, nrow)),
                "psCDR3.treg" = unlist(lapply(treg.psCDR3, nrow)),
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
write.csv(subset.stats, paste(fig.file.path, "subset.stats.csv", sep = "/"))

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
heatmap.2(mat_data, scale = "row", tracecol=NA, dendrogram = "none", Rowv = FALSE, Colv = FALSE,
          col = my_palette)
dev.off()

# Notice interesting patters. Many neighborhoods restricted to a single compartment.  
# Need to get sequences that contribute to this pattern

teff.spec.nbrs <- neighborhood.df$neighborhood[neighborhood.df$samp.count.teff > 0 & neighborhood.df$samp.count.treg == 0]
treg.spec.nbrs <- neighborhood.df$neighborhood[neighborhood.df$samp.count.treg > 0 & neighborhood.df$samp.count.teff == 0]

# Get the sequences
# Look for CDR3s of top neighborhoods in Teff
teff.spec.seq <- lapply(teff.parse, function(x){
  CDR3s <- unlist(core.neighborhoods[as.integer(teff.spec.nbrs)])
  # Iterate through all of the neighborhoods
  seqs <- x$CDR3.amino.acid.sequence[x$CDR3.amino.acid.sequence %in% CDR3s]
  return(seqs)
}) %>% unlist()

treg.spec.seq <- lapply(treg.parse, function(x){
  CDR3s <- unlist(core.neighborhoods[as.integer(treg.spec.nbrs)])
  # Iterate through all of the neighborhoods
  seqs <- x$CDR3.amino.acid.sequence[x$CDR3.amino.acid.sequence %in% CDR3s]
  return(seqs)
}) %>%  unlist()

# Make a graph with the Teff/Treg information
pairs.graph <- graph_from_data_frame(pairs)
# Add in node attributes based on Teff or Treg presence
pairs.graph <- igraph::as_data_frame(pairs.graph, what = "both")
pairs.graph$vertices$subset <- apply(pairs.graph$vertices, 1, function(x){
  subset <- ifelse(x[["name"]] %in% teff.spec.seq, "teff", ifelse(x[["name"]] %in% treg.spec.seq, "treg", "none"))
})

# Want to make less binary, determine if any are "skewed" to a subset
pairs.graph$vertices$teff.to.treg <- unlist(lapply(pairs.graph$vertices$name, function(x){
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
  prop <- teff.prop.mean / treg.prop.mean
  if(is.nan(prop)){
    return(0)
  } else if(is.infinite(prop)) {
    return(50)
  } else{
    return(prop)
  }

}))

# Want to make less binary, determine if any are "skewed" to a subset
pairs.graph$vertices$treg.to.teff <- unlist(lapply(pairs.graph$vertices$name, function(x){
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
  prop <- treg.prop.mean / teff.prop.mean
  if(is.nan(prop)){
    return(0)
  } else if(is.infinite(prop)) {
    return(50)
  } else{
    return(prop)
  }
  
}))

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


