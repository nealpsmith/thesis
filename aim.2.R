rm(list = ls())
library(plyr) # load before dplyr to avoid conflicts due to masking
library(dplyr)
library(tcR)
library(ggplot2)
library(grid)
library(ggrepel)
library(magrittr)
library(stringi)
library(stringr)
library(foreach)
library(doParallel)
library(data.table)
source("tcR_functions.R") #includes several modifications to tcR


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

enriched.CDR3s.df <- do.call(rbind, enriched.CDR3s)
### Motif analysis ###

# Create list of 4mers and 5mers
threemer_list <- nmer.gen(enriched.CDR3s.df$CDR3.amino.acid.sequence, 3)
fourmer_list <- nmer.gen(enriched.CDR3s.df$CDR3.amino.acid.sequence, 4)
fivemer_list <- nmer.gen(enriched.CDR3s.df$CDR3.amino.acid.sequence, 5)

# Determine rates of nmers in enriched clones
threemer.enriched <- nmer.search(threemer_list, enriched.CDR3s.df, count.col = 2)
fourmer.enriched <- nmer.search(fourmer_list, enriched.CDR3s.df, count.col = 2)
fivemer.enriched <- nmer.search(fivemer_list, enriched.CDR3s.df, count.col = 2)

# Get the negative data
CD154.neg.CDR3s <- lapply(data.parse[grep("neg", names(data.parse))], function(x){
  x = select(x, CDR3.amino.acid.sequence, Read.count)
  return(x)
}) %>% do.call(rbind, .)

# Aggregate all repeated CDR3s, adding up read counts
CD154.neg.CDR3s <- aggregate(. ~ CDR3.amino.acid.sequence, data = CD154.neg.CDR3s, sum)

# Set cores for parallel processing
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

# Use paralell processing to look for the nmers in the CD154- CDR3s
threemer.neg <- rbindlist(foreach(i = 1:length(threemer.enriched$nmer)) %dopar% {
  source(paste(raw.file.path, "neals_tcr_functions.R", sep = "/"))
  data = nmer.search(threemer.enriched$nmer[i], CD154.neg.CDR3s, CDR3.col = 1, count.col = 2)
  data
})

# Set cores for parallel processing
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

fourmer.neg <- rbindlist(foreach(i = 1:length(fourmer.enriched$nmer)) %dopar% {
  source(paste(raw.file.path, "neals_tcr_functions.R", sep = "/"))
  data = nmer.search(fourmer.enriched$nmer[i], CD154.neg.CDR3s, CDR3.col = 1, count.col = 2)
  data
})

# Set cores for parallel processing
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

fivemer.neg <- rbindlist(foreach(i = 1:length(fivemer.enriched$nmer)) %dopar% {
  source(paste(raw.file.path, "neals_tcr_functions.R", sep = "/"))
  data = nmer.search(fivemer.enriched$nmer[i], CD154.neg.CDR3s, CDR3.col = 1, count.col = 2)
  data
})

stopCluster(cl)

# Need to get some statistics for negative data
data.stats <- cloneset.stats(data.parse)

# Put the data into lists
motif.pos <- list(threemer = threemer.enriched,
                  fourmer = fourmer.enriched,
                  fivemer = fivemer.enriched)
motif.neg <- list(threemer = threemer.neg,
                  fourmer = fourmer.neg,
                  fivemer = fivemer.neg)


for(i in 1:length(motif.pos)){
  motif.pos[[i]]$rate.neg <- motif.neg[[i]]$rate
  motif.pos[[i]]$enrichment <- motif.pos[[i]]$rate / motif.neg[[i]]$rate
  df <- data.frame(enrich.count = motif.pos[[i]]$clone.count,
                   neg.count = motif.neg[[i]]$clone.count,
                   enrich.total = sum(enriched.CDR3s.df$activated.count),
                   neg.total = sum(data.stats[,"Sum.reads"][grep("neg", rownames(data.stats))]))
  motif.pos[[i]]$p.value <- apply(df, 1, function(x){
    gt <- GTest(matrix(x, ncol = 2))
    return(gt$p.value)
  })
}

motif.pos <- lapply(motif.pos, function(x){
  x$enrichment[x$enrichment == "Inf" | x$enrichment > 500] <- 500
  return(x)
})

# Do FDR adjustment
motif.pos <- lapply(motif.pos, function(x) FDR.function(x, 0.05))

pdf(paste(fig.file.path, "nmer.plots.pdf", sep = "/"), 12, 7)
ggplot(data = motif.pos$fourmer, aes(x = rate.neg, y = rate)) +
  geom_jitter(data = motif.pos$fourmer[!motif.pos$fourmer$FDRsig == 1 |
                                         !motif.pos$fourmer$rate > 0.005 |
                                         !motif.pos$fourmer$rate > 2 * motif.pos$fourmer$rate.neg &
                                         motif.pos$fourmer$unique.CDR3s > 2,], color = "grey", size = 1) +
  geom_jitter(data = motif.pos$fourmer[motif.pos$fourmer$FDRsig == 1 &
                                         motif.pos$fourmer$rate > 0.005 &
                                         motif.pos$fourmer$rate > 2 * motif.pos$fourmer$rate.neg &
                                         motif.pos$fourmer$unique.CDR3s > 2,],
              color = "red", aes(size = unique.CDR3s)) +
  geom_text_repel(data = motif.pos$fourmer[motif.pos$fourmer$FDRsig == 1 &
                                             motif.pos$fourmer$rate >0.005 &
                                             motif.pos$fourmer$rate > 2 * motif.pos$fourmer$rate.neg &
                                             motif.pos$fourmer$unique.CDR3s > 2,],
                  aes(label = nmer),
                  min.segment.length = 0.12, size = 5) +
  scale_x_continuous(limits = c(0.000, 0.0125)) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Proportion in CD154- CDR3s") + ylab("Proportion in enriched CDR3s") +
  labs(size = "Unique CDR3s") +
  ggtitle("4mers in Enriched CDR3s vs. CD154- CDR3s") +
  theme_bw(base_size = 18)


ggplot(data = motif.pos$fivemer, aes(x = rate.neg, y = rate)) +
  geom_jitter(data = motif.pos$fivemer[!motif.pos$fivemer$FDRsig == 1 |
                                         !motif.pos$fivemer$rate > motif.pos$fivemer$rate.neg |
                                         !motif.pos$fivemer$unique.CDR3s > 2 |
                                         !motif.pos$fivemer$rate >=0.0025,], color = "grey", size = 1) +
  geom_jitter(data = motif.pos$fivemer[motif.pos$fivemer$FDRsig == 1 &
                                         motif.pos$fivemer$rate > motif.pos$fivemer$rate.neg &
                                         motif.pos$fivemer$unique.CDR3s > 2 &
                                         motif.pos$fivemer$rate >= 0.0025,],
              color = "red", aes(size = unique.CDR3s)) +
  geom_text_repel(data = motif.pos$fivemer[motif.pos$fivemer$FDRsig == 1 &
                                             motif.pos$fivemer$rate > motif.pos$fivemer$rate.neg &
                                             motif.pos$fivemer$unique.CDR3s > 2 &
                                             motif.pos$fivemer$rate > 0.0025,],
                  aes(label = nmer),
                  min.segment.length = 0.05) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Proportion in CD154- CDR3s") + ylab("Proportion in enriched CDR3s") +
  labs(size = "Unique CDR3s") +
  ggtitle("5mers in Enriched CDR3s vs. CD154- CDR3s") +
  theme_bw(base_size = 18)
dev.off()

# Want to look at how many people contribute to each nmer in their enriched
motif.pos <- lapply(motif.pos, function(x){
  x$unique.people <- apply(x, 1, function(y){
    CDR3.df <- enriched.CDR3s.df[grepl(y["nmer"], substr(enriched.CDR3s.df$CDR3.amino.acid.sequence, 4, 
                                                         nchar(enriched.CDR3s.df$CDR3.amino.acid.sequence) - 3)),]
    people <- length(unique(CDR3.df$id))
    return(people)
  })
  return(x)
})

# Look at the distributions of unique CDR3s, enrichments, and people who contribute to an nmer
# First unique CDR3s
pdf(paste(fig.file.path, "unique.nmer.CDR3s.hist.pdf", sep = "/"), 10, 7)
ggplot(data=motif.pos$threemer[motif.pos$threemer$unique.people > 1,], aes(unique.CDR3s)) + geom_histogram(bins = 100) +
  xlab("Unique CDR3s") + ylab("nmer count") +
  ggtitle("Unique CDR3s distribution: 3mers") +
  theme_bw()
ggplot(data=motif.pos$fourmer[motif.pos$fourmer$unique.people > 1,], aes(unique.CDR3s)) + geom_histogram(bins = 100) +
  xlab("Unique CDR3s") + ylab("nmer count") +
  ggtitle("Unique CDR3s distribution: 4mers") +
  theme_bw()
ggplot(data=motif.pos$fivemer[motif.pos$fivemer$unique.people > 1,], aes(unique.CDR3s)) + geom_histogram(bins = 54) +
  xlab("Unique CDR3s") + ylab("nmer count") +
  ggtitle("Unique CDR3s distribution: 5mers") +
  theme_bw()
dev.off()

# Next look at enrichment (proportion of nmer in enriched / proportion of nmer in CD154-)
pdf(paste(fig.file.path, "enrichment.nmers.hist.pdf", sep = "/"), 10, 7)
ggplot(data = motif.pos$threemer[motif.pos$threemer$unique.people > 1,], aes(enrichment)) + geom_histogram(bins = 100) +
  xlab("prop.enriched / prop.resting") + ylab("nmer count") +
  ggtitle("enrichment distribution: 3mers") +
  theme_bw()
ggplot(data = motif.pos$fourmer[motif.pos$fourmer$unique.people > 1,], aes(enrichment)) + geom_histogram(bins = 100) +
xlab("prop.enriched / prop.resting") + ylab("nmer count") +
  ggtitle("enrichment distribution: 4mers") +
  theme_bw()
ggplot(data = motif.pos$fivemer[motif.pos$fivemer$unique.people > 1,], aes(enrichment)) + geom_histogram(bins = 100) +
xlab("prop.enriched / prop.resting") + ylab("nmer count") +
  ggtitle("enrichment distribution: 5mers") +
  theme_bw()
dev.off()

pdf(paste(fig.file.path, "unique.people.nmer.hist.pdf", sep = "/"), 10, 7)
ggplot(data = motif.pos$threemer, aes(unique.people)) + geom_histogram(bins = 20) +
  xlab("Unique people") + ylab("nmer count") +
  ggtitle("unique people distribution: 3mers") +
  theme_bw()
ggplot(data = motif.pos$fourmer, aes(unique.people)) + geom_histogram(bins = 17) +
  xlab("Unique people") + ylab("nmer count") +
  ggtitle("unique people distribution: 4mers") +
  theme_bw()
ggplot(data = motif.pos$fivemer, aes(unique.people)) + geom_histogram(bins = 17) +
  xlab("Unique people") + ylab("nmer count") +
  ggtitle("unique people distribution: 5mers") +
  theme_bw()
dev.off()


# Make plots looking at enrichment vs. unique people
pdf(paste(fig.file.path, "nmer.plots.by.enrichment.pdf", sep = "/"), 10, 7)
ggplot(data = motif.pos$fourmer, aes(x = unique.people, y = log(enrichment))) + 
  geom_jitter(data = motif.pos$fourmer[!motif.pos$fourmer$unique.CDR3s >= 3 |
                                         ! motif.pos$fourmer$unique.people >=3 |
                                         ! motif.pos$fourmer$enrichment >= 10,], color = "grey", size = 1) +
  geom_jitter(data = motif.pos$fourmer[motif.pos$fourmer$unique.CDR3s >= 3 &
                                         motif.pos$fourmer$unique.people >=3 &
                                         motif.pos$fourmer$enrichment >= 10,], color = "red", aes(size = unique.CDR3s)) +
  scale_size_continuous(range = c(2,6)) +
  geom_text_repel(data = motif.pos$fourmer[motif.pos$fourmer$unique.CDR3s >= 3 &
                                             motif.pos$fourmer$unique.people >=3 &
                                             motif.pos$fourmer$enrichment >= 10,], aes(label = nmer)) +
  xlab("Unique subjects") + ylab("Log(enrichment)") +
  theme_bw()

ggplot(data = motif.pos$fivemer, aes(x = unique.people, y = log(enrichment))) + 
  geom_jitter(data = motif.pos$fivemer[!motif.pos$fivemer$unique.CDR3s >= 3 |
                                         ! motif.pos$fivemer$unique.people >=3 |
                                         ! motif.pos$fivemer$enrichment >= 10,], color = "grey", size = 1) +
  geom_jitter(data = motif.pos$fivemer[motif.pos$fivemer$unique.CDR3s >= 3 &
                                         motif.pos$fivemer$unique.people >=3 &
                                         motif.pos$fivemer$enrichment >= 10,], color = "red", aes(size = unique.CDR3s)) +
  scale_size_continuous(range = c(2,6)) +
  geom_text_repel(data = motif.pos$fivemer[motif.pos$fivemer$unique.CDR3s >= 3 &
                                             motif.pos$fivemer$unique.people >=3 &
                                             motif.pos$fivemer$enrichment >= 10,], aes(label = nmer)) +
  xlab("Unique subjects") + ylab("Log(enrichment)") +
  theme_bw()
dev.off()

# Look at motifs in all CD154+ vs. Cd154-

# Get all CD154+ data
CD154.pos.CDR3s <- lapply(data.parse[grep("pos", names(data.parse))], function(x){
  x = select(x, CDR3.amino.acid.sequence, Read.count)
  return(x)
}) %>% do.call(rbind, .)


# Look for nmers in all CD154+ data
threemer.pos.all <- nmer.search(threemer.enriched$nmer, CD154.pos.CDR3s)
fourmer.pos.all <- nmer.search(fourmer.enriched$nmer, CD154.pos.CDR3s)
fivemer.pos.all <- nmer.search(fivemer.enriched$nmer, CD154.pos.CDR3s)

# Put the data into a list
motif.pos.all <- list(threemer = threemer.pos.all,
                  fourmer = fourmer.pos.all,
                  fivemer = fivemer.pos.all)

# Add this in with the other information
for(i in 1:length(motif.pos)){
  motif.pos[[i]]$rate.pos.all <- motif.pos.all[[i]]$rate
}


pdf(paste(fig.file.path, "enriched.nmers.all.cd154.pos.vs.neg.pdf", sep = "/"), 10, 7)
ggplot(data = motif.pos$fourmer, aes(x = rate.neg, y = rate.pos.all)) +
  geom_jitter(data = motif.pos$fourmer[!motif.pos$fourmer$FDRsig == 1 |
                                         !motif.pos$fourmer$rate > 0.005 |
                                         !motif.pos$fourmer$rate > 2 * motif.pos$fourmer$rate.neg &
                                         motif.pos$fourmer$unique.CDR3s > 2,], color = "grey", size = 1) +
  geom_jitter(data = motif.pos$fourmer[motif.pos$fourmer$FDRsig == 1 &
                                         motif.pos$fourmer$rate > 0.005 &
                                         motif.pos$fourmer$rate > 2 * motif.pos$fourmer$rate.neg &
                                         motif.pos$fourmer$unique.CDR3s > 2,],
              color = "red", aes(size = unique.CDR3s)) +
  geom_text_repel(data = motif.pos$fourmer[motif.pos$fourmer$FDRsig == 1 &
                                             motif.pos$fourmer$rate >0.005 &
                                             motif.pos$fourmer$rate > 2 * motif.pos$fourmer$rate.neg &
                                             motif.pos$fourmer$unique.CDR3s > 2,],
                  aes(label = nmer),
                  min.segment.length = 0.12, size = 5) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Proportion in CD154- CDR3s") + ylab("Proportion in enriched CDR3s") +
  labs(size = "Unique CDR3s") +
  ggtitle("4mers in Enriched CDR3s vs. CD154- CDR3s") +
  theme_bw(base_size = 18)


ggplot(data = motif.pos$fivemer, aes(x = rate.neg, y = rate.pos.all)) +
  geom_jitter(data = motif.pos$fivemer[!motif.pos$fivemer$FDRsig == 1 |
                                         !motif.pos$fivemer$rate > motif.pos$fivemer$rate.neg |
                                         !motif.pos$fivemer$unique.CDR3s > 2 |
                                         !motif.pos$fivemer$rate >=0.0025,], color = "grey", size = 1) +
  geom_jitter(data = motif.pos$fivemer[motif.pos$fivemer$FDRsig == 1 &
                                         motif.pos$fivemer$rate > motif.pos$fivemer$rate.neg &
                                         motif.pos$fivemer$unique.CDR3s > 2 &
                                         motif.pos$fivemer$rate >= 0.0025,],
              color = "red", aes(size = unique.CDR3s)) +
  geom_text_repel(data = motif.pos$fivemer[motif.pos$fivemer$FDRsig == 1 &
                                             motif.pos$fivemer$rate > motif.pos$fivemer$rate.neg &
                                             motif.pos$fivemer$unique.CDR3s > 2 &
                                             motif.pos$fivemer$rate > 0.0025,],
                  aes(label = nmer),
                  min.segment.length = 0.05) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Proportion in CD154- CDR3s") + ylab("Proportion in enriched CDR3s") +
  labs(size = "Unique CDR3s") +
  ggtitle("5mers in Enriched CDR3s vs. CD154- CDR3s") +
  theme_bw(base_size = 18)

dev.off()


# Based on these distributions, make cutoffs and see what percentile they represent
top.threemer <- motif.pos$threemer[motif.pos$threemer$unique.CDR3s >= 3 & 
                          motif.pos$threemer$enrichment >= 10 &
                          motif.pos$threemer$unique.people >=3,]

# Determine percentile of total 3mers this is
nrow(top.threemer) / nrow(motif.pos$threemer)

# Based on these distributions, make cutoffs and see what percentile they represent
top.fourmer <- motif.pos$fourmer[motif.pos$fourmer$unique.CDR3s >= 3 & 
                                     motif.pos$fourmer$enrichment >= 10 &
                                     motif.pos$fourmer$unique.people >=3,]
# Determine percentile
nrow(top.fourmer) / nrow(motif.pos$fourmer)


# Based on these distributions, make cutoffs and see what percentile they represent
top.fivemer <- motif.pos$fivemer[motif.pos$fivemer$unique.CDR3s >= 3 & 
                                   motif.pos$fivemer$enrichment >= 10 &
                                   motif.pos$fivemer$unique.people >=3,]
# Determine percentile
nrow(top.fivemer) / nrow(motif.pos$fivemer)


### Look at nmers in a set of activated CDR3s from the infectious disease realm ### 
infec.data <- read.table("infectious.disease.tsv", sep = "\t", header = TRUE)
infec.data$CDR3 <- as.character(infec.data$CDR3)
infec.data$count <- 1

# Look for nmers in infectious disease data
threemer.infec <- nmer.search(threemer.enriched$nmer, infec.data, CDR3.col = 3, 
                                count.col = ncol(infec.data))
fourmer.infec <- nmer.search(fourmer.enriched$nmer, infec.data, CDR3.col = 3, 
                              count.col = ncol(infec.data))
fivemer.infec <- nmer.search(fivemer.enriched$nmer, infec.data, CDR3.col = 3, 
                              count.col = ncol(infec.data))

motif.pos.infec <- list(threemer = threemer.enriched,
                        fourmer = fourmer.enriched,
                        fivemer = fivemer.enriched)
motif.infec.list <- list(threemer = threemer.infec,
                         fourmer = fourmer.infec,
                         fivemer = fivemer.infec)

for(i in 1:length(motif.pos.infec)){
  motif.pos.infec[[i]]$rate.infec <- motif.infec.list[[i]]$rate
}

# Need to adjust the "rate" for enriched CDR3s
# Don't have read counts for infectious disease so we will not look at counts!  Will treat every CDR3 the same

for(i in 1:length(motif.pos.infec)){
  motif.pos.infec[[i]]$rate <- motif.pos.infec[[i]]$unique.CDR3s / 
    length(unique(enriched.CDR3s.df$CDR3.amino.acid.sequence))
}

pdf(paste(fig.file.path, "nmer.plots.enriched.vs.infectious.pdf", sep = "/"), 10, 7 )
ggplot(data = motif.pos.infec$fourmer, aes(x = rate.infec, y = rate)) +
  geom_jitter(data = motif.pos.infec$fourmer[!motif.pos.infec$fourmer$rate > 0.005 |
                                         !motif.pos.infec$fourmer$rate > 2 * motif.pos.infec$fourmer$rate.infec &
                                         motif.pos.infec$fourmer$unique.CDR3s > 2,], color = "grey", size = 1) +
  geom_jitter(data = motif.pos.infec$fourmer[motif.pos.infec$fourmer$rate > 0.005 &
                                         motif.pos.infec$fourmer$rate > 2 * motif.pos.infec$fourmer$rate.infec &
                                         motif.pos.infec$fourmer$unique.CDR3s > 2,],
              color = "red", aes(size = unique.CDR3s)) +
  geom_text_repel(data = motif.pos.infec$fourmer[motif.pos.infec$fourmer$rate > 0.005 &
                                                   motif.pos.infec$fourmer$rate > 2 * motif.pos.infec$fourmer$rate.infec &
                                                   motif.pos.infec$fourmer$unique.CDR3s > 2,],
                  aes(label = nmer),
                   size = 5) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Proportion in infectious CDR3s") + ylab("Proportion in enriched CDR3s") +
  scale_x_continuous(limits = c(0, 0.02)) +
  labs(size = "Unique CDR3s") +
  ggtitle("4mers ininfectious disease CDR3s vs. peanut enriched CDR3s") +
  theme_bw(base_size = 18)

ggplot(data = motif.pos.infec$fivemer, aes(x = rate.infec, y = rate)) +
  geom_jitter(data = motif.pos.infec$fivemer[!motif.pos.infec$fivemer$rate > 0.005 |
                                               !motif.pos.infec$fivemer$rate > 2 * motif.pos.infec$fivemer$rate.infec &
                                               motif.pos.infec$fivemer$unique.CDR3s > 2,], color = "grey", size = 1) +
  geom_jitter(data = motif.pos.infec$fivemer[motif.pos.infec$fivemer$rate > 0.005 &
                                               motif.pos.infec$fivemer$rate > 2 * motif.pos.infec$fivemer$rate.infec &
                                               motif.pos.infec$fivemer$unique.CDR3s > 2,],
              color = "red", aes(size = unique.CDR3s)) +
  geom_text_repel(data = motif.pos.infec$fivemer[motif.pos.infec$fivemer$rate > 0.005 &
                                                   motif.pos.infec$fivemer$rate > 2 * motif.pos.infec$fivemer$rate.infec &
                                                   motif.pos.infec$fivemer$unique.CDR3s > 2,],
                  aes(label = nmer),
                  min.segment.length = 0.12, size = 5) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Proportion in infectious CDR3s") + ylab("Proportion in enriched CDR3s") +
  scale_x_continuous(limits = c(0, 0.015)) +
  labs(size = "Unique CDR3s") +
  ggtitle("5mers ininfectious disease CDR3s vs. peanut enriched CDR3s") +
  theme_bw(base_size = 18)
dev.off()
