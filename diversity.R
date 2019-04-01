library(alakazam)
library(dplyr)
library(magrittr)
library(ggplot2)
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
       raw.file.path <- "~/data/thesis/neals.thesis.data"},
       stop("I don't recognize your username, type Sys.info() to find out what it is.")
)



# Load in the data
load(paste(raw.file.path, "parsed.data.baseline.rda", sep = "/"))
load(paste(raw.file.path, "enriched.CDR3s.rda", sep = "/"))

enriched.CDR3s.df <- do.call(rbind, enriched.CDR3s)

seqs <- apply(enriched.CDR3s.df, 1, function(x){
  unlist(rep(x["CDR3.amino.acid.sequence"], as.numeric(x["activated.count"])), use.names = FALSE)}) %>%
  unlist(., use.names = FALSE)

ids <- apply(enriched.CDR3s.df, 1, function(x){
  unlist(rep(x["id"], as.numeric(x["activated.count"])), use.names = FALSE)}) %>%
  unlist(., use.names = FALSE)

diversity.df <- data.frame(CDR3.amino.acid.sequence = seqs,
                           group = ids)
diversity.df <- diversity.df[diversity.df$group == "100",]

# Set up comparison within each subject
ids <- unique(enriched.CDR3s.df$id)

# Make a dataframe for each subject ID
all.data.by.id <- lapply(ids, function(x){
  # Get the enriched df
  enriched.df<- enriched.CDR3s[[x]] %>%
    select(., CDR3.amino.acid.sequence, activated.count) %>%
    `colnames<-`(c("CDR3.amino.acid.sequence", "Read.count"))
  enriched.df$group <- "psCDR3"
  # Get all CD154+ data
  all.pos.df <- data.parse[[paste(x, "pos", sep = "_")]] %>%
    select(., CDR3.amino.acid.sequence, Read.count)
  all.pos.df$group <- "all_CD154+"
  # Get all the CD154- data
  neg.df <- data.parse[[paste(x, "neg", sep = "_")]] %>%
    select(., CDR3.amino.acid.sequence, Read.count)
  neg.df$group <- "CD154-"
  
  all.df <- do.call(rbind, list(enriched.df, all.pos.df, neg.df))
  return(all.df)
  
}) %>% `names<-`(ids)

# Calculate diversity for each subject
diversity <- lapply(all.data.by.id, function(x){
  div <- rarefyDiversity(x, group = "group", clone = "CDR3.amino.acid.sequence",
                         copy = "Read.count")
})

div.sig <- lapply(all.data.by.id, function(x){
  sig.tests <- list()
  # Iterate all values of q
  for( i in 0:2){
    test <- testDiversity(x, group = "group", clone = "CDR3.amino.acid.sequence",
                          copy = "Read.count", q = i)
    p.vals <- test@tests
    p.vals$q <- i
    sig.tests[[as.character(i)]] <- p.vals
  }
  sig.tests <- do.call(rbind, sig.tests)
  
  return(sig.tests)
}) %>% do.call(rbind, .)



pdf(paste(fig.file.path ,"diversity.curves.psCDR3s.pdf", sep = "/"))
for(i in 1:length(diversity)){
  plot(diversity[[i]], 
       legend_title="Sample", main_title = paste("diversity:", names(diversity[i])))
  
}
dev.off()

# Look at fold-change in diversity from CD154-:CD154+, CD154-:psCDR3, CD154+:psCDR3
div.fc <- lapply(all.data.by.id, function(x){
  fold.changes <- list()
  # Iterate all values of q
  for( i in 0:2){
    test <- testDiversity(x, group = "group", clone = "CDR3.amino.acid.sequence",
                          copy = "Read.count", q = i)
    df <- data.frame(comp = c("psCDR3toCD154neg", "psCDR3toCD154pos", "CD154negtoCD154pos"),
                     fold.change = c(test@summary["CD154-",]$MEAN / test@summary["psCDR3",]$MEAN,
                                     test@summary["all_CD154+",]$MEAN / test@summary["psCDR3",]$MEAN,
                                     test@summary["CD154-",]$MEAN / test@summary["all_CD154+",]$MEAN))
    df$q <- i
    fold.changes[[as.character(i)]] <- df
  }
  fold.changes <- do.call(rbind, fold.changes)
  
  return(fold.changes)
}) %>% do.call(rbind, .)

pdf(paste(fig.file.path, "fold.change.at.div.indeces.pdf", sep = "/"), 16, 12)
ggplot(div.fc, aes(x = comp, y = fold.change)) + geom_boxplot() +
  scale_y_continuous(limits = c(0, 8)) +
  facet_wrap(~q, scales = "free_y") +
  ylab("Fold change") + xlab("") +
  scale_x_discrete(labels = c("CD154-:CD154+", "CD154-:psCDR3", "CD154+:psCDR3")) +
  ggtitle("fold change of diversity at different indeces") +
  theme_bw(base_size = 15)
dev.off()