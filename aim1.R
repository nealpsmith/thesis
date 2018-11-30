library(plyr) # load before dplyr to avoid conflicts due to masking
library(dplyr)
library(tcR)
library(ggplot2)
library(grid)
library(gridBase)
library(ggrepel)
library(magrittr)
library(stringi)
library(stringr)
library(RVAideMemoire)
library(DescTools)
library(entropy) # Masks entropy function in tcR
source("tcR_functions.R") #includes several modifications to tcR


# Get file path for output figures and raw data files
fig.file.path <- "C:/Users/nealp/Dropbox (Personal)/Extension School/Thesis/figures"
raw.file.path <- "C:/Users/nealp/Dropbox (Partners HealthCare)/Projects/PNOIT2-1037/TCRB sequencing and HLA typing data"

source(paste(raw.file.path, "neals_tcr_functions.R", sep = "/"))

# Add clinical phenotype data
# Includes 4 hyporeactives who were confirmed allergic by challenge!
reactive.ids <- c("69", "80", "84", "85", "89", "93", "100", "104", "106", "107", "27",
                  "81", "101", "94")

# Collect all of the baseline data
orig.data.path <- "C:/Users/nealp/Dropbox (Partners HealthCare)/PAID-UP, PNOIT2 Reactive vs Non-reactive paper/TCR Sequencing"

# Get original TCRseq data and limit to just confirmed allergic
orig.data.pos <- list.files(paste(orig.data.path, "raw_data_CD154+", sep = "/"))[
  grep(paste(reactive.ids, sep = "", collapse = "|"), 
       list.files(paste(orig.data.path, "raw_data_CD154+", sep = "/")))]

orig.data.neg <- list.files(paste(orig.data.path, "raw_data_CD154-CD69-", sep = "/"))[
  grep(paste(reactive.ids, sep = "", collapse = "|"), 
       list.files(paste(orig.data.path, "raw_data_CD154+", sep = "/")))]

# Get newer files from in-house runs
new.files.plate.1 <- list.files(paste(raw.file.path, "plate.1", sep = "/"))[
  grep("BL", list.files(paste(raw.file.path, "plate.1", sep = "/")))
]
new.files.plate.2 <- list.files(paste(raw.file.path, "plate.2", sep = "/"))[
  grep("BL", list.files(paste(raw.file.path, "plate.2", sep = "/")))
]
new.files.plate.3 <- list.files(paste(raw.file.path, "plate.3", sep = "/"))[
  grep("BL", list.files(paste(raw.file.path, "plate.3", sep = "/")))
  ]

# Parse the data 
data.parse <- parse.adaptive.file.list(c(paste(orig.data.path, "raw_data_CD154+", orig.data.pos, sep = "/"), 
                                         paste(orig.data.path, "raw_data_CD154-CD69-", orig.data.neg, sep = "/"),
                                         paste(raw.file.path, "plate.1", new.files.plate.1, sep = "/"),
                                         paste(raw.file.path, "plate.2", new.files.plate.2, sep = "/"),
                                         paste(raw.file.path, "plate.3", new.files.plate.3, sep = "/")))

# Adjust the names of the files to make them consistent
names(data.parse)[grep("\\d[0-9][B]", names(data.parse))] <- gsub("B_TCRB", "_neg",
                                                                  names(data.parse[grep("\\d[0-9][B]",
                                                                                        names(data.parse))]))

names(data.parse)[grep("_TCRB", names(data.parse))] <- gsub("_TCRB", "_pos", names(data.parse[grep("_TCRB",
                                                                                                   names(data.parse))]))

names(data.parse)[grep("_BL_", names(data.parse))] <- gsub("_BL_", "_", names(data.parse[grep("_BL_",
                                                                                              names(data.parse))]))
save(data.parse, file = paste(raw.file.path, "neals.thesis.data/parsed.data.baseline.rda", sep = "/"))

# Look for enriched CDR3s
enriched.CDR3s <- CDR3.enriched(data.parse, pos.pattern = "_pos", neg.pattern = "_neg")
save(enriched.CDR3s, file = paste(raw.file.path, "neals.thesis.data/enriched.CDR3s.rda", sep = "/"))

enriched.CDR3.df <- do.call(rbind, enriched.CDR3s)

# Hamming analysis
# Perform hamming on entire list of peanut CDR3s
min.ham <- CDR3_ham(enriched.CDR3.df$CDR3.amino.acid.sequence)


# Create a list of all of the CDR3s from CD154- to be used as a control
rest.CDR3s <- lapply(data.parse[grep("neg",names(data.parse))], function(x){
  x <- unique(x$CDR3.amino.acid.sequence)
}) %>% unlist(., use.names = FALSE)

# Get hamming distance on all CD154+ and CD154-
neg.ham.df <- data.frame(distance = 0:20)

for(i in 1:100){
  neg.ham <- CDR3_ham(sample(rest.CDR3s, nrow(enriched.CDR3.df)))
  neg.ham.df <- left_join(neg.ham.df, neg.ham, by = "distance", all = TRUE)
}

# Clean up negative dataframe, add median percent and standard deviation
neg.ham.df[is.na(neg.ham.df)] <- 0
neg.ham.df$median <- apply(neg.ham.df[2:ncol(neg.ham.df)], 1, median)
neg.ham.df$sd <- apply(neg.ham.df[2:ncol(neg.ham.df)], 1, sd)

# Create a list of all of the CDR3s from CD154- to be used as a control
act.CDR3s <- lapply(data.parse[grep("pos",names(data.parse))], function(x){
  x <- unique(x$CDR3.amino.acid.sequence)
}) %>% unlist(., use.names = FALSE)


# Do the same for all CD154+
pos.all.ham.df <- data.frame(distance = 0:20)
for(i in 1:100){
  pos.all.ham <- CDR3_ham(sample(act.CDR3s, nrow(enriched.CDR3.df)))
  pos.all.ham.df <- left_join(pos.all.ham.df, pos.all.ham, by = "distance", all = TRUE) 
}

# Clean up negative dataframe, add median percent and standard deviation
pos.all.ham.df[is.na(pos.all.ham.df)] <- 0
pos.all.ham.df$median <- apply(pos.all.ham.df[2:ncol(pos.all.ham.df)], 1, median)
pos.all.ham.df$sd <- apply(pos.all.ham.df[2:ncol(pos.all.ham.df)], 1, sd)


sd <-  neg.ham.df$sd[as.numeric(as.character(neg.ham.df$distance)) <=15] %>%
  append(pos.all.ham.df$sd[as.numeric(as.character(pos.all.ham.df$distance)) <=15]) %>%
  append(rep(NA, 16))

min.ham <- full_join(min.ham, neg.ham.df, by = "distance") %>%
  select(c(1, 2, ncol(.) - 1)) %>%
  full_join(., pos.all.ham.df, by = "distance") %>%
  select(c(1,3, ncol(.) - 1, 2))
colnames(min.ham) <- c("distance", "CD154-", "CD154+", "enriched.CD154+")

min.ham[is.na(min.ham)] <- 0
min.ham <- min.ham[!min.ham$distance == "Inf",]
min.ham <- min.ham[order(min.ham$distance),]

min.ham.melt <- melt(min.ham, id.vars = "distance")
min.ham.melt$distance <- as.integer(min.ham.melt$distance)
min.ham.melt <- min.ham.melt[min.ham.melt$distance <=15,] %>%
  cbind(., sd)

pdf(paste(fig.file.path, "hamming.distance.peanut.CDR3s.pdf", sep = "/"), 10, 7)
ggplot(min.ham.melt, aes(x = distance, y = value, fill = variable, width = .80)) +
  geom_bar(stat = "identity", position = position_dodge()) + scale_x_continuous(breaks = c(0:15)) +
  scale_fill_manual(values = c('#E69F00', "#999999", "#b30000"), labels = c("CD154- T cells", "All CD154+ T cells",
                                                                            "Peanut Enriched CD154+ cells")) +
  geom_errorbar(aes(ymin = value, ymax = sd + value), width = .3, position = position_dodge(0.8)) +
  ylab("percent of total clones") + xlab("AA differences to next closest CDR3") +
  ggtitle("Hamming distance for peanut CDR3s")
dev.off()


# shannon entropy
# # Test: prove to yourself the above function is actually calculating the Shannon entropy
# num <- sum(ant.bind$activated.value)
# x <- apply(ant.bind, 1, function(x){
#   freq = as.numeric(x["activated.value"]) / num
#   val = freq * log2(freq)
#   return(val)
# })
# # This should equal the same as the calculated value above
# -sum(x)
# 

# Calculate the entropy of random selection of CDR3s
# First need to create a vector of all CDR3s, repeating sequences based on counts, so they can be sampled randomly
all.CDR3.reads <- lapply(data.parse, function(x){
  CDR3s <- apply(x, 1, function(y){
    seq <- rep(y["CDR3.amino.acid.sequence"], as.numeric(y["Read.count"]))
    return(seq)
  }) %>% unlist(., use.names = FALSE)
}) %>% unlist(., use.names = FALSE)

# # Test: Make sure length of the vector is the same as the sum of all reads
# sum(unlist(lapply(baseline.combined, function(x) sum(x$Read.count)))) == length(all.CDR3.reads)

# Shannon index by individual
# Compare the shannon index of the antigen binding region from a single individual and compare that to their own CD154-
enrich.split.bind <- lapply(enriched.CDR3s, function(x){
  x$CDR3.amino.acid.sequence <- tcr_trim(x$CDR3.amino.acid.sequence, 3, 3)
  x <- aggregate(activated.count ~ CDR3.amino.acid.sequence, data = x, FUN = sum)
  return(x)
})


enrich.split.bind.entr <- lapply(enrich.split.bind, function(x) entropy(x$activated.count, unit = "log2"))



# Look at the CD154- for each individual
# Need to create random vectors of CDR3s for each individual
CD154.neg.CDR3.reads <- lapply(data.parse[grep("neg", names(data.parse))], function(x){
  CDR3s <- apply(x, 1, function(y){
    seq <- rep(y["CDR3.amino.acid.sequence"], as.numeric(y["Read.count"]))
    return(seq)
  }) %>% unlist(., use.names = FALSE)
})
# Change the names on the list of CD154- reads to be just subject IDs
names(CD154.neg.CDR3.reads) <- substr(names(CD154.neg.CDR3.reads), 0, nchar(names(CD154.neg.CDR3.reads)) - 4)

# Create a blank dataframe to store entropy data
neg.entropy.df <- data.frame(id = names(CD154.neg.CDR3.reads),
                             median.entropy.neg = numeric(length = length(CD154.neg.CDR3.reads)),
                             sd.neg = numeric(length = length(CD154.neg.CDR3.reads)))

# Iterate through the negative reads from each individual
for(i in 1:length(CD154.neg.CDR3.reads)){
  # Create empty vector
  entr.vec <- c()
  # Resample negative reads 50 times
  for(j in 1:50){
    # Get corresponding enriched data sample
    pos.sample <- as.data.frame(enrich.split.bind[[match(names(CD154.neg.CDR3.reads[i]), names(enrich.split.bind))]])
    # Sample negative data at the same amount as the corresponding enriched counts
    sample <- sample(CD154.neg.CDR3.reads[[i]], sum(pos.sample$activated.count)) %>%
      # Trim to likely-antigen contact region
      tcr_trim(., 3, 3) %>%
      # Summarize the data
      table() %>% as.data.frame()
    # Determine entropy of the sampling
    entropy <- entropy(sample$Freq, unit = "log2")
    # Add it to a vector
    entr.vec <- c(entr.vec, entropy)
  }
  # Add median of 50 resamples to dataframe
  neg.entropy.df[i,c(2, 3)] <- c(median(entr.vec), sd(entr.vec))
  
}

# Include total CD154+ in entropy measurments
# Need to create random vectors of CDR3s for each individual of their CD154+
CD154.pos.CDR3.reads <- lapply(data.parse[grep("pos", names(data.parse))], function(x){
  CDR3s <- apply(x, 1, function(y){
    seq <- rep(y["CDR3.amino.acid.sequence"], as.numeric(y["Read.count"]))
    return(seq)
  }) %>% unlist(., use.names = FALSE)
})

# # TEST: Make sure there are the correct number of CDR3 sequences in the vectors
# length(unlist(CD154.pos.CDR3.reads)) == sum(unlist(lapply(data.parse[grep("pos", names(data.parse))],
#                                                           function(x){
#                                                             sum(x$Read.count)
#                                                             })))

# Change names to be just IDs
names(CD154.pos.CDR3.reads) <- substr(names(CD154.pos.CDR3.reads), 0, nchar(names(CD154.pos.CDR3.reads)) - 4)

# Create a blank dataframe to store entropy data
pos.entropy.df <- data.frame(id = names(CD154.pos.CDR3.reads),
                             median.entropy.pos = numeric(length = length(CD154.pos.CDR3.reads)),
                             sd.pos = numeric(length = length(CD154.pos.CDR3.reads)))

# Iterate through all of the CD154 positive data from each subject
for(i in 1:length(CD154.pos.CDR3.reads)){
  # Create empty vector for resampling
  entr.vec <- c()
  # Resample CD154+ data 50 times
  for(j in 1:50){
    # Get corresponding enriched data
    pos.sample <- as.data.frame(enrich.split.bind[[match(names(CD154.pos.CDR3.reads[i]), names(enrich.split.bind))]])
    # Sample the CD154+ CDR3s in the same amount as the corresponding enriched reads
    sample <- sample(CD154.pos.CDR3.reads[[i]], sum(pos.sample$activated.count)) %>%
      # Trim CDR3s to likely-antigen contact region
      tcr_trim(., 3, 3) %>%
      # Summarize the data
      table() %>% as.data.frame()
    # Determine entropy
    entropy <- entropy(sample$Freq, unit = "log2")
    # Add it to the vector
    entr.vec <- c(entr.vec, entropy)
  }
  # Put median entropy into the dataframe
  pos.entropy.df[i,c(2, 3)] <- c(median(entr.vec), sd(entr.vec))
  
}

# Create dataframe for plotting
plot.df <- data.frame(id = names(enrich.split.bind.entr),
                      enrich.entropy = unlist(enrich.split.bind.entr)) %>%
  left_join(., neg.entropy.df[c(1,2)], by = "id") %>%
  left_join(., pos.entropy.df[c(1,2)], by = "id") %>%
  select(., c(1,3,4,2)) %>%
  melt(id.var = "id")


# Plot the data
pdf(paste(fig.file.path, "shannons.index.with.CD154.pos.pdf", sep = "/"), 7, 7)
ggplot(plot.df, aes(x = variable, y = value, group = id)) + geom_line(size = 1) +
  geom_point(size = 3) +
  ylab("Shannon's H index (bits)") + xlab("") +
  scale_x_discrete(labels = c("CD154- CDR3s","CD154+ CDR3s", "Enriched CDR3s")) +
  ggtitle("Shannon's index by subject") +
  theme_bw(base_size = 14)
dev.off()


# Downsampling of enriched CDR3s to compare clonality of enriched populations between individuals
# Create vectors of CDR3s
enrich.vector.list <- lapply(enriched.CDR3s, function(x){
  CDR3s <- apply(x, 1, function(y){
    seq <- rep(y["CDR3.amino.acid.sequence"], as.numeric(y["activated.count"]))
    return(seq)
  }) %>% unlist(., use.names = FALSE)
  CDR3s <- tcr_trim(CDR3s, 3, 3)
})

# Determine sample sizes to resample with
sample.sizes <- seq(10, 100, 10)

# create a dataframe to store resample info in
downsample.df <- data.frame(id = names(enrich.vector.list),
                 matrix(ncol = length(sample.sizes), nrow = length(enrich.vector.list)))

# Add column names to the dataframe
colnames(downsample.df) <- c("id", sample.sizes)

# Iterate through all the sample sizes
for(i in sample.sizes){
  # Determine number of resamplings to do
  num.resamples <- 100
  # Create a dataframe to store info after every iteration
  resample.df <- data.frame(id = names(enrich.vector.list),
                            matrix(ncol = num.resamples, nrow = length(enrich.vector.list)))
  # Perform resampling
  for(j in 1:num.resamples){
    # Determine entropy for every subject
    samp.entropy <- lapply(enrich.vector.list, function(x){
      # If sample size larger than size of subject's library, ignore it
      if(length(x) < i){
        entropy <- NA
      } else{
        # Get sample and limit to antigen-contact region
        sample <- sample(x, i) %>%
          table() %>% as.data.frame()
        # Determine entropy of sample
        entropy <- entropy(sample$Freq, unit = "log2")
      }
    })
    # Add info about entropy for every iteration
    resample.df[,j + 1] <- unlist(samp.entropy)
  }
  # Determine the median for all resamplings
  resample.df$median <- apply(resample.df[2:ncol(resample.df)], 1, median)
  # Add medians to resample dataframe
  downsample.df[,as.character(i)] <- resample.df$median
}

# Make a dataframe for plotting
plot.df <- melt(downsample.df, id.vars = "id")

# Plot the data
pdf(paste(fig.file.path, "shannons.index.enriched.downsampling.pdf", sep = "/"), 10, 7)
ggplot(plot.df, aes(x = variable, y = value, group = id)) + geom_line() + geom_point() +
  xlab("sample size") + ylab("Shannon H index (bits)") +
  ggtitle("Downsampling of enriched CDR3s: Shannon's entropy") +
  theme_bw(base_size = 14)
dev.off()
