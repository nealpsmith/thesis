# library(plyr) # load before dplyr to avoid conflicts due to masking
library(dplyr)
# library(tcR)
# library(ggplot2)
# library(grid)
# library(ggrepel)
library(magrittr)
# library(stringi)
# library(stringr)
# library(data.table)
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

# Load in the raw data
# Load in the data
load(paste(raw.file.path, "parsed.data.baseline.rda", sep = "/"))
load(paste(raw.file.path, "enriched.CDR3s.rda", sep = "/"))

# Make enriched data into single dataframe
enriched.CDR3s.df <- do.call(rbind, enriched.CDR3s)

# Get the negative data
CD154.neg.CDR3s <- lapply(data.parse[grep("neg", names(data.parse))], function(x){
  x = select(x, CDR3.amino.acid.sequence, Read.count)
  return(x)
}) %>% do.call(rbind, .)

# Look for discontinuous motifs
disc.fourmers.enriched <- disc.4mers(enriched.CDR3s.df, CDR3.col = "CDR3.amino.acid.sequence",
                                      count.col = "activated.count")

# Add in read proportion to dataframe
disc.fourmers.enriched$prop = disc.fourmers.enriched$count / sum(enriched.CDR3s.df$activated.count)

### For CD154- CDR3s, do parallel processing to make things go faster ###
# Set cores for parallel processing
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

# Iterate across all of the CDR3s
disc.fourmers.neg <- foreach(i = 1:nrow(CD154.neg.CDR3s), .combine = rbind) %dopar% {
  source("neals_tcr_functions.R")
  # Make sure CDR3 is long enough
  if(nchar(CD154.neg.CDR3s[i,]$CDR3.amino.acid.sequence) >=10){
    # look at discontinuous 4mer on an individual CDR3 basis
    data <- disc.4mers(CD154.neg.CDR3s[i,], CDR3.col = "CDR3.amino.acid.sequence",
                       count.col = "Read.count")
    data
  }
  
}
stopCluster(cl)
remove(cores); remove(cl)

disc.fourmers.neg <- aggregate(. ~ nmer, data = disc.fourmers.neg, sum)

disc.fourmers.neg$prop <- disc.fourmers.neg$count / sum(CD154.neg.CDR3s$Read.count)

# Join negative data to positive data
disc.fourmers.enriched <- left_join(disc.fourmers.enriched, disc.fourmers.neg, by = "nmer")
colnames(disc.fourmers.enriched) <- c("nmer", "unique.CDR3s.enriched", "activated.count",
                                       "prop.enriched", "unique.CDR3s.neg", "resting.count",
                                       "prop.neg")
disc.fourmers.enriched$enrichment <- disc.fourmers.enriched$prop.enriched / disc.fourmers.enriched$prop.neg

pdf("discontinous.4mers.pdf", 12, 7)
ggplot(data = disc.fourmers.enriched, aes(x = prop.neg, y = prop.enriched)) +
  geom_jitter(data = disc.fourmers.enriched[!disc.fourmers.enriched$enrichment > 2 &
                                         disc.fourmers.enriched$unique.CDR3s.enriched >= 1,], color = "grey", size = 1) +
  geom_jitter(data = disc.fourmers.enriched[disc.fourmers.enriched$enrichment >=2 &
                                               disc.fourmers.enriched$unique.CDR3s.enriched > 2 &
                                               disc.fourmers.enriched$prop.enriched >= 0.005,],
              color = "red", aes(size = unique.CDR3s.enriched)) +
  geom_text_repel(data = disc.fourmers.enriched[disc.fourmers.enriched$enrichment >=2 &
                                                   disc.fourmers.enriched$unique.CDR3s.enriched > 2 &
                                                   disc.fourmers.enriched$prop.enriched >=0.005,],
                  aes(label = nmer),
                  min.segment.length = 0.12, size = 5) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Proportion in CD154- CDR3s") + ylab("Proportion in enriched CDR3s") +
  labs(size = "Unique CDR3s") +
  ggtitle("discontinuous 4mers in Enriched CDR3s vs. CD154- CDR3s") +
  theme_bw(base_size = 18)
dev.off()


# Set cores for parallel processing
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

# Iterate across all of the CDR3s
disc.fourmers.enriched$unique.people <- unlist(foreach(i = 1:length(disc.fourmers.enriched$nmer))
                                                %dopar% {
          source(neals.func.path)
          people <- unique.people.discontinuous(disc.fourmers.enriched$nmer[i],
                                                enriched.CDR3s.df, "CDR3.amino.acid.sequence",
                                                "id", motif.size = 4)
          people
      }
)


top.disc.fourmers <- disc.fourmers.enriched[disc.fourmers.enriched$unique.CDR3s.enriched >=3 &
                               disc.fourmers.enriched$enrichment >=10 &
                               disc.fourmers.enriched$unique.people >=3,]
save(file = paste(raw.file.path, "top.disc.fourmers.rda", sep = "/"), top.disc.fourmers)


### Look at discontinuous 5mers ###
# Look for discontinuous motifs
disc.fivemers.enriched <- disc.5mers(enriched.CDR3s.df, CDR3.col = "CDR3.amino.acid.sequence",
                                      count.col = "activated.count")
# Add in read proportion to dataframe
disc.fivemers.enriched$prop = disc.fivemers.enriched$count / sum(enriched.CDR3s.df$activated.count)

### For CD154- CDR3s, do parallel processing to make things go faster ###
# Set cores for parallel processing
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

# Iterate across all of the CDR3s
disc.fivemers.neg <- foreach(i = 1:nrow(CD154.neg.CDR3s), .combine = rbind) %dopar% {
  source("neals_tcr_functions.R")
  # Make sure CDR3 is long enough
  if(nchar(CD154.neg.CDR3s[i,]$CDR3.amino.acid.sequence) >=11){
    # look at discontinuous 4mer on an individual CDR3 basis
    data <- disc.5mers(CD154.neg.CDR3s[i,], CDR3.col = "CDR3.amino.acid.sequence",
                       count.col = "Read.count")
    data
  }
  
}
stopCluster(cl)
remove(cores); remove(cl)



disc.fivemers.neg <- aggregate(. ~ nmer, data = disc.fivemers.neg, sum)

disc.fivemers.neg$prop <- disc.fivemers.neg$count / sum(CD154.neg.CDR3s$Read.count)

# Join negative data to positive data
disc.fivemers.enriched <- left_join(disc.fivemers.enriched, disc.fivemers.neg, by = "nmer")
colnames(disc.fivemers.enriched) <- c("nmer", "unique.CDR3s.enriched", "activated.count",
                                      "prop.enriched", "unique.CDR3s.neg", "resting.count",
                                      "prop.neg")
disc.fivemers.enriched$enrichment <- disc.fivemers.enriched$prop.enriched / disc.fivemers.enriched$prop.neg


# Set cores for parallel processing
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

# Iterate across all of the CDR3s
disc.fivemers.enriched$unique.people <- unlist(foreach(i = 1:length(disc.fivemers.enriched$nmer))
                                               %dopar% {
                                                 source(paste(raw.file.path, "neals_tcr_functions.R", sep = "/"))
                                                 people <- unique.people.discontinuous(disc.fivemers.enriched$nmer[i],
                                                                                       enriched.CDR3s.df, "CDR3.amino.acid.sequence",
                                                                                       "id")
                                                 people
                                               }
)

top.disc.fivemers <- disc.fivemers.enriched[disc.fivemers.enriched$unique.CDR3s.enriched >=3 &
                                              disc.fivemers.enriched$enrichment >=10 &
                                              disc.fivemers.enriched$unique.people >=3,]
save(file = paste(raw.file.path, "top.disc.fivemers.rda", sep = "/"), top.disc.fivemers)

