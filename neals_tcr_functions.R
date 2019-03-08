# Useful functions created by Neal for TCRseq analysis

# Function to get middle of TCR trimming by any parameters you want.
# Clone_vec is character vector of CDR3s to be trimmed.
# x is distance to trim from front, y is distance to trim from back
tcr_trim <- function(clone_vec, x, y){
    trimmed_CDR3s <- c()
    for(clones in 1:length(clone_vec)){
        middle <- stri_sub(clone_vec[[clones]], x, -y)
        trimmed_CDR3s <- c(trimmed_CDR3s, middle)
    }
    trimmed_CDR3s;
}


# Efficient function to generate nmers in the antigen-binding region of CDR3s of any size
# CDR3 = character vector of CDR3s, len = length of nmer
nmer.gen <- function(CDR3, len){
    all.nmers <- sapply(as.character(CDR3), function(x){
        mid <- substr(x, 4, nchar(x) - 3)
        nmers <- substring(mid, seq(1, nchar(mid) - (len - 1)), seq(len, nchar(mid)))
        return(nmers)

    })
    nmers <- unique(unlist(all.nmers, use.names = FALSE))
    nmers <- nmers[nchar(nmers) == len]

    return(nmers)
}


# Function to determine number of CDR3s that have a specific nmer
# nmer is list of nmer, CDR3s is list of CDR3s
unique_CDR3s <- function(nmer, trimmed_CDR3s){
    count_list <- c()
    for(i in 1:length(nmer)){
        x <- length(unique(grep(nmer[[i]], trimmed_CDR3s, value = TRUE)))
        count_list <- c(count_list, x)
    }
    data.frame(nmer = nmer,
                unique.CDR3s = count_list);

}

# A flexible function that takes in a list of nmers and will return the number of counts, the number of unique CDR3s and
# the number of clones have that nmer.  Will return a dataframe.
nmer.search  <- function(nmer, df, CDR3.col = 1, count.col = 2){
    df.out <- data.frame(nmer = nmer,
                     tot.clones = numeric(length = length(nmer)),
                     unique.CDR3s = numeric(length = length(nmer)),
                     clone.count = numeric(length = length(nmer)),
                     rate = numeric(length = length(nmer)))

    list <- lapply(nmer, function(nmer){
        df.nmer <- df[grepl(nmer, substr(df[[CDR3.col]], 4, nchar(df[[CDR3.col]]) - 3)) == TRUE,]
        info <- list(length(df.nmer[[CDR3.col]]),length(unique(df.nmer[[CDR3.col]])), sum(df.nmer[[count.col]]),
                     sum(df.nmer[[count.col]]) / sum(df[[count.col]]))
        return(info)

    })
    df.out$tot.clones <- sapply(list, '[[', 1)
    df.out$unique.CDR3s <- sapply(list, '[[', 2)
    df.out$clone.count <- sapply(list, '[[', 3)
    df.out$rate <- sapply(list, '[[', 4)
    df.out;
}

# Function to get the r-squared value for a given set of data (2mer, 3mer ect.) df = nmer summary chart w/pos and neg data
lm_eqn <- function(df){
    m <- lm( stim.count ~ median.count.CD154neg, df);
    eq <- substitute(italic(r)^2~"="~r2,
                     list(r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
}

# Calculate chi-squared values
chi_square <- function(df, x, y){
    temp <- data.frame(stim.hits = df[[x]],
                       median.count.CD154.neg = df[[y]])
    df1 <- apply(temp, 1, function(x) {
        ch <- chisq.test(matrix(x, ncol =2, byrow = TRUE))
        ch$p.value})
    df <- cbind(df, df1);

}

# Function to calculate hamming distance of a vector of CDR3 sequences
# CDR3s = character vector of CDR3s, method is either hamming or levenshtein
CDR3_dist <- function(CDR3s, method){
    library(stringdist)
    ham.table <- lapply(seq_along(CDR3s), function(x){
        min(stringdist(CDR3s[x], CDR3s[-x], method = method))
    }) %>% unlist(., use.names = FALSE) %>% table(.) %>%
        as.data.frame(.) %>%
        transform(., percent = Freq / length(CDR3s) * 100)
    ham.table$Freq <- NULL
    colnames(ham.table) <- c("distance", "percent")
    ham.table$distance <- as.numeric(as.character(ham.table$distance))
    return(ham.table)
}

# Create a function to use an FDR cutoff to evaluate enrichment
FDR.function <- function(df, cutoff){
    df$q.value <- p.adjust(df$p.value, method = "BH")
    df$FDRsig <- as.numeric(df$q.value <= cutoff)
    return(df)
}

# Function to analyze CD154+/- data and return enriched CDR3 dataframe
CDR3.enriched <- function(data, pos.pattern, neg.pattern, FDR.cutoff = 0.05){
    # Make a list of pairs of dataframes for each sample
    pair.list <- vector("list", length = length(data) / 2)

    # Get the data files in pairs: CD154+ and CD154-
    for(i in seq(1, length(data), 1)){
        pair <- c(names(data[i]),
                  names(data)[gsub(neg.pattern, "", names(data)[i]) ==
                                               gsub(pos.pattern, "", names(data))])
        pair.list[[i]] <- pair
    }
    pair.list <- pair.list[lengths(pair.list) == 2]

    # Put all sample sets in same order for consistency
    pair.list <- lapply(pair.list, function(x) x <- x[c(grep(pos.pattern, x), grep(neg.pattern, x))])

    # Put each pair of data through statistical enrichment
    parse.data <- lapply(pair.list, function(x){
        # Get Statistics (for total counts)
        stats <- cloneset.stats(data[x])
        # Determine counts of each CDR3 found in both CD154+ and CD154-
        parse <- shared.repertoire(data[x], .type = 'a0rc', .min.ppl = 1, .verbose = F)
        parse$People <- NULL
        # Adjust the column names to make the next code easier and to ensure data can be "rbind"ed afterwards
        colnames(parse) <- c("CDR3.amino.acid.sequence", "activated.count", "resting.count")
        # Get rid of any CDR3s that do not show up in the CD154+ data
        parse <- parse[!is.na(parse$activated.count),]
        # Make the count for any that did not show up in the CD154- 0
        parse[is.na(parse$resting.count),]$resting.count <- 0
        # Make columns for total CD154+ and CD154- reads
        parse$activated.reads <- stats[,8][grep(pos.pattern, rownames(stats))]
        parse$resting.reads <- stats[,8][grep(neg.pattern, rownames(stats))]
        # Add an ID column based on the leading numeric characters in the names of the files (should be subject/sample ID)
        parse$id <- sub("\\D*(\\d+).*", "\\1", x[1])
        parse$timepoint <- substr(x[1], 0,nchar(x[1]) - nchar(pos.pattern)) %>%
            sub("^[0-9]+[_]", "", .)


        # Determine p.values based on G test
        # Make a new dataframe for g-test
        g.test.df <- parse[2:5]

        # Adjust read counts for contingency table
        g.test.df$activated.reads <- g.test.df$activated.reads - g.test.df$activated.count
        g.test.df$resting.reads <- g.test.df$resting.reads - g.test.df$resting.count

        # Determine p-values for every CDR3
        parse$p.value <- apply(g.test.df, 1, function(x){
            gt <- GTest(matrix(x, ncol = 2))
            return(gt$p.value)
        })

        # Adjust p-values with FDR cutoff determined by user
        parse <- FDR.function(parse, FDR.cutoff) %>%
            # Filter out those that do not meet FDR cutoff or count criteria
            filter(FDRsig == 1 & activated.count > 1 &
                       ((activated.count / activated.reads) /
                            (resting.count / resting.reads)) > 1)

        return(parse)
    })
    names <- sapply(pair.list, "[[", 1)
    names(parse.data) <-substr(names, 0, nchar(names) - nchar(pos.pattern))
    return(parse.data)

}

# Function to loook for discontinuous 4mers
# df = dataframe with CDR3s and counts, CDR3.col = index or character vector for the CDR3 column
# count.col = index or character vector for count column
disc.4mers <- function(df, CDR3.col, count.col){
  # Exclude CDR3s that are too short
  CDR3.df <- df[nchar(df[[CDR3.col]]) >= 10,]

  # Iterate through the CDR3s
  motifs <- apply(CDR3.df, 1, function(x){
    # Trim the CDR3
    trimmed.CDR3 <- substr(x[[CDR3.col]], 4, nchar(x[[CDR3.col]]) - 3)
    # Make an empty vector
    combos <- c()
    # Loop across a CDR3 and get all possible discontinuous 3mers
    for(i in 1:(nchar(trimmed.CDR3) - 3)){
      comb.1 <- paste(substr(trimmed.CDR3, i, i), "X", substr(trimmed.CDR3, i + 2, i + 3), sep = "")
      comb.2 <- paste(substr(trimmed.CDR3, i, i + 1), "X", substr(trimmed.CDR3, i + 3, i + 3), sep = "")
      combos <- c(combos, comb.1, comb.2)
    }
    # Only count a discontinuous motif once per CDR3
    combos <- unique(combos)
    # Get the counts of the CDR3 where these nmers were derived
    counts <- as.integer(x[[count.col]])
    return(data.frame("nmer" = combos, "unique.CDR3s" = 1, "count" = counts))
  })
   #combo.df <- data.frame("nmer" = unlist(motifs, use.names = FALSE), "unique.CDR3s" = 1)
  combo.df <- do.call(rbind, motifs)
  combo.df <- aggregate(. ~ nmer, data = combo.df, sum)
  return(combo.df)
}

disc.5mers <- function(df, CDR3.col, count.col){
  # Exclude CDR3s that are too short
  CDR3.df <- df[nchar(df[[CDR3.col]]) >= 11,]
  
  # Iterate through the CDR3s
  motifs <- apply(CDR3.df, 1, function(x){
    # Trim the CDR3
    trimmed.CDR3 <- substr(x[[CDR3.col]], 4, nchar(x[[CDR3.col]]) - 3)
    # Make an empty vector
    combos <- c()
    # Loop across a CDR3 and get all possible discontinuous 3mers
    for(i in 1:(nchar(trimmed.CDR3) - 4)){
      comb.1 <- paste(substr(trimmed.CDR3, i, i), "X", substr(trimmed.CDR3, i + 2, i + 4), sep = "")
      comb.2 <- paste(substr(trimmed.CDR3, i, i + 1), "X", substr(trimmed.CDR3, i + 3, i + 4), sep = "")
      comb.3 <- paste(substr(trimmed.CDR3, i, i + 2), "X", substr(trimmed.CDR3, i + 4, i + 4), sep = "")
      comb.4 <- paste(substr(trimmed.CDR3, i, i), "XX", substr(trimmed.CDR3, i + 3, i + 4), sep = "")
      comb.5 <- paste(substr(trimmed.CDR3, i, i + 1), "XX", substr(trimmed.CDR3, i + 4, i + 4), sep = "")
      comb.6 <- paste(substr(trimmed.CDR3, i, i), "X", substr(trimmed.CDR3, i + 2, i + 2), "X",
                      substr(trimmed.CDR3, i + 4, i + 4), sep = "")
      combos <- c(combos, comb.1, comb.2, comb.3, comb.4, comb.5, comb.6)
    }
    # Only count a discontinuous motif once per CDR3
    combos <- unique(combos)
    # Get the counts of the CDR3 where these nmers were derived
    counts <- as.integer(x[[count.col]])
    return(data.frame("nmer" = combos, "unique.CDR3s" = 1, "count" = counts))
  })
  #combo.df <- data.frame("nmer" = unlist(motifs, use.names = FALSE), "unique.CDR3s" = 1)
  combo.df <- do.call(rbind, motifs)
  combo.df <- aggregate(. ~ nmer, data = combo.df, sum)
  return(combo.df)
}

unique.people.discontinuous <- function(nmer, df, CDR3.col, id.col, motif.size){
  if(motif.size == 4){
    # Get rid of CDR3s that are too short to contribute
    CDR3.df <- df[nchar(df[[CDR3.col]]) >= 10,]
    
    # Iterate across the CDR3s looking for which ones have an nmer
    ids <- apply(CDR3.df, 1, function(x){
      # Trim the CDR3
      trimmed.CDR3 <- substr(x[[CDR3.col]], 4, nchar(x[[CDR3.col]]) - 3)
      # Make an empty vector
      combos <- c()
      # Loop across a CDR3 and get all possible discontinuous 3mers
      for(i in 1:(nchar(trimmed.CDR3) - 3)){
        comb.1 <- paste(substr(trimmed.CDR3, i, i), "X", substr(trimmed.CDR3, i + 2, i + 3), sep = "")
        comb.2 <- paste(substr(trimmed.CDR3, i, i + 1), "X", substr(trimmed.CDR3, i + 3, i + 3), sep = "")
        combos <- c(combos, comb.1, comb.2)
      }
      # Only count a discontinuous motif once per CDR3
      combos <- unique(combos)
      if(nmer %in% combos){
        return(x[[id.col]])
      }
    })
    ids <- Filter(Negate(is.null), ids)
    return(length(unique(ids)))
  }
  if(motif.size == 5){
    # Exclude CDR3s that are too short
    CDR3.df <- df[nchar(df[[CDR3.col]]) >= 11,]
    # Iterate across the CDR3s looking for which ones have an nmer
    ids <- apply(CDR3.df, 1, function(x){
      # Trim the CDR3
      trimmed.CDR3 <- substr(x[[CDR3.col]], 4, nchar(x[[CDR3.col]]) - 3)
      # Make an empty vector
      combos <- c()
      # Loop across a CDR3 and get all possible discontinuous 3mers
      for(i in 1:(nchar(trimmed.CDR3) - 4)){
        comb.1 <- paste(substr(trimmed.CDR3, i, i), "X", substr(trimmed.CDR3, i + 2, i + 4), sep = "")
        comb.2 <- paste(substr(trimmed.CDR3, i, i + 1), "X", substr(trimmed.CDR3, i + 3, i + 4), sep = "")
        comb.3 <- paste(substr(trimmed.CDR3, i, i + 2), "X", substr(trimmed.CDR3, i + 4, i + 4), sep = "")
        comb.4 <- paste(substr(trimmed.CDR3, i, i), "XX", substr(trimmed.CDR3, i + 3, i + 4), sep = "")
        comb.5 <- paste(substr(trimmed.CDR3, i, i + 1), "XX", substr(trimmed.CDR3, i + 4, i + 4), sep = "")
        comb.6 <- paste(substr(trimmed.CDR3, i, i), "X", substr(trimmed.CDR3, i + 2, i + 2), "X",
                        substr(trimmed.CDR3, i + 4, i + 4), sep = "")
        combos <- c(combos, comb.1, comb.2, comb.3, comb.4, comb.5, comb.6)
      }
      # Only count a discontinuous motif once per CDR3
      combos <- unique(combos)
      if(nmer %in% combos){
        return(x[[id.col]])
      }
    })
    ids <- Filter(Negate(is.null), ids)
    return(length(unique(ids)))
  }
  
}
# Function that creates edge between two CDR3s that share a discontinuous motif
find_pairs_disc <- function(nmer, df, CDR3.col, motif.size){
  if(motif.size == 4){
    # Get rid of CDR3s that are too short to contribute
    CDR3.df <- df[nchar(df[[CDR3.col]]) >= 10,]
    
    # Iterate across the CDR3s looking for which ones have an nmer
    CDR3s <- apply(CDR3.df, 1, function(x){
      # Trim the CDR3
      trimmed.CDR3 <- substr(x[[CDR3.col]], 4, nchar(x[[CDR3.col]]) - 3)
      # Make an empty vector
      combos <- c()
      # Loop across a CDR3 and get all possible discontinuous 3mers
      for(i in 1:(nchar(trimmed.CDR3) - 3)){
        comb.1 <- paste(substr(trimmed.CDR3, i, i), "X", substr(trimmed.CDR3, i + 2, i + 3), sep = "")
        comb.2 <- paste(substr(trimmed.CDR3, i, i + 1), "X", substr(trimmed.CDR3, i + 3, i + 3), sep = "")
        combos <- c(combos, comb.1, comb.2)
      }
      # Only count a discontinuous motif once per CDR3
      combos <- unique(combos)
      if(nmer %in% combos){
        return(x[[CDR3.col]])
      }
    })
    CDR3s <- Filter(Negate(is.null), CDR3s)
    
    CDR3s <- unlist(CDR3s, use.names = FALSE)
    if(length(CDR3s) >= 2){
      pairs <- combn(CDR3s, 2)
      df <- data.frame("from.cdr3" = pairs[1,],
                       "to.cdr3" = pairs[2,])
      return(df)
    } else{
      return(data.frame("from.cdr3" = NULL,
                        "to.cdr3" = NULL))
    }
  }
  
  if(motif.size == 5){
    # Get rid of CDR3s that are too short to contribute
    CDR3.df <- df[nchar(df[[CDR3.col]]) >= 11,]
    
    CDR3s <- apply(CDR3.df, 1, function(x){
      # Trim the CDR3
      trimmed.CDR3 <- substr(x[[CDR3.col]], 4, nchar(x[[CDR3.col]]) - 3)
      # Make an empty vector
      combos <- c()
      # Loop across a CDR3 and get all possible discontinuous 3mers
      for(i in 1:(nchar(trimmed.CDR3) - 4)){
        comb.1 <- paste(substr(trimmed.CDR3, i, i), "X", substr(trimmed.CDR3, i + 2, i + 4), sep = "")
        comb.2 <- paste(substr(trimmed.CDR3, i, i + 1), "X", substr(trimmed.CDR3, i + 3, i + 4), sep = "")
        comb.3 <- paste(substr(trimmed.CDR3, i, i + 2), "X", substr(trimmed.CDR3, i + 4, i + 4), sep = "")
        comb.4 <- paste(substr(trimmed.CDR3, i, i), "XX", substr(trimmed.CDR3, i + 3, i + 4), sep = "")
        comb.5 <- paste(substr(trimmed.CDR3, i, i + 1), "XX", substr(trimmed.CDR3, i + 4, i + 4), sep = "")
        comb.6 <- paste(substr(trimmed.CDR3, i, i), "X", substr(trimmed.CDR3, i + 2, i + 2), "X",
                        substr(trimmed.CDR3, i + 4, i + 4), sep = "")
        combos <- c(combos, comb.1, comb.2, comb.3, comb.4, comb.5, comb.6)
      }
      # Only count a discontinuous motif once per CDR3
      combos <- unique(combos)
      if(nmer %in% combos){
        return(x[[CDR3.col]])
      }
    })
    CDR3s <- Filter(Negate(is.null), CDR3s)
    
    CDR3s <- unlist(CDR3s, use.names = FALSE)
    if(length(CDR3s) >= 2){
      pairs <- combn(CDR3s, 2)
      df <- data.frame("from.cdr3" = pairs[1,],
                       "to.cdr3" = pairs[2,])
      return(df)
    } else{
      return(data.frame("from.cdr3" = NULL,
                        "to.cdr3" = NULL))
    }
    
  }
  
}

# Get pairs based on continuous nmers
find_pairs_cont <- function(nmer, df, CDR3.col){
  # Iterate across the CDR3s looking for which ones have an nmer
  CDR3s <- df[[CDR3.col]][grepl(nmer, substr(df[[CDR3.col]], 4, nchar(df[[CDR3.col]]) - 3))]
  if(length(CDR3s) >= 2){
    pairs <- combn(CDR3s, 2)
    df <- data.frame("from.cdr3" = pairs[1,],
                     "to.cdr3" = pairs[2,])
    return(df)
  } else{
    return(data.frame("from.cdr3" = NULL,
                      "to.cdr3" = NULL))
  }
}

# Get pairs based on continuous nmers
find_cont <- function(nmer, df, CDR3.col){
  # Iterate across the CDR3s looking for which ones have an nmer
  CDR3s <- df[[CDR3.col]][grepl(nmer, substr(df[[CDR3.col]], 4, nchar(df[[CDR3.col]]) - 3))]
  return(CDR3s)
}


# Function that will return CDR3s with a particular discontinuous nmer
find_disc <- function(nmer, CDR3s){
  if(nchar(nmer) == 4){
    # Get rid of CDR3s that are too short to contribute
    CDR3s <- CDR3s[nchar(CDR3s) >= 10]
    # Iterate across the CDR3s looking for which ones have an nmer
    CDR3s <- sapply(CDR3s, function(x){
      # Trim the CDR3
      trimmed.CDR3 <- substr(x, 4, nchar(x) - 3)
      # Make an empty vector
      combos <- c()
      # Loop across a CDR3 and get all possible discontinuous 3mers
      for(i in 1:(nchar(trimmed.CDR3) - 3)){
        comb.1 <- paste(substr(trimmed.CDR3, i, i), "X", substr(trimmed.CDR3, i + 2, i + 3), sep = "")
        comb.2 <- paste(substr(trimmed.CDR3, i, i + 1), "X", substr(trimmed.CDR3, i + 3, i + 3), sep = "")
        combos <- c(combos, comb.1, comb.2)
      }
      # Only count a discontinuous motif once per CDR3
      combos <- unique(combos)
      if(nmer %in% combos){
        return(x)
      }
    })
    CDR3s <- as.character(Filter(Negate(is.null), CDR3s))
    
    return(CDR3s)
  }
  if(nchar(nmer) == 5){
    # Get rid of CDR3s that are too short to contribute
    CDR3s <- CDR3s[nchar(CDR3s) >= 11]
    
    # Iterate across the CDR3s looking for which ones have an nmer
    CDR3s <- sapply(CDR3s, function(x){
      # Trim the CDR3
      trimmed.CDR3 <- substr(x, 4, nchar(x) - 3)
      # Make an empty vector
      combos <- c()
      # Loop across a CDR3 and get all possible discontinuous 3mers
      for(i in 1:(nchar(trimmed.CDR3) - 3)){
        comb.1 <- paste(substr(trimmed.CDR3, i, i), "X", substr(trimmed.CDR3, i + 2, i + 4), sep = "")
        comb.2 <- paste(substr(trimmed.CDR3, i, i + 1), "X", substr(trimmed.CDR3, i + 3, i + 4), sep = "")
        comb.3 <- paste(substr(trimmed.CDR3, i, i + 2), "X", substr(trimmed.CDR3, i + 4, i + 4), sep = "")
        comb.4 <- paste(substr(trimmed.CDR3, i, i), "XX", substr(trimmed.CDR3, i + 3, i + 4), sep = "")
        comb.5 <- paste(substr(trimmed.CDR3, i, i + 1), "XX", substr(trimmed.CDR3, i + 4, i + 4), sep = "")
        comb.6 <- paste(substr(trimmed.CDR3, i, i), "X", substr(trimmed.CDR3, i + 2, i + 2), "X",
                        substr(trimmed.CDR3, i + 4, i + 4), sep = "")
        combos <- c(combos, comb.1, comb.2, comb.3, comb.4, comb.5, comb.6)
      }
      # Only count a discontinuous motif once per CDR3
      combos <- unique(combos)
      if(nmer %in% combos){
        return(x)
      }
    })
    CDR3s <- as.character(Filter(Negate(is.null), CDR3s))
    return(CDR3s)
    
  }
}

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



disc.fourmer.gen  <- function(df, CDR3.col){
  # Get rid of CDR3s that are too short to contribute
  CDR3.df <- df[nchar(df[[CDR3.col]]) >= 10,]
  
  # create all discontinuous nmers for the dataframe
  discs <- sapply(CDR3.df[[CDR3.col]], function(x){
    # Trim the CDR3
    trimmed.CDR3 <- substr(x, 4, nchar(x) - 3)
    
    # Make an empty vector
    combos <- c()
    
    # Loop across a CDR3 and get all possible discontinuous 4mers
    for(j in 1:(nchar(trimmed.CDR3) - 3)){
      comb.1 <- paste(substr(trimmed.CDR3, j, j), "X", substr(trimmed.CDR3, j + 2, j + 3), sep = "")
      comb.2 <- paste(substr(trimmed.CDR3, j, j + 1), "X", substr(trimmed.CDR3, j + 3, j + 3), sep = "")
      combos <- c(combos, comb.1, comb.2)
    }
    # Only count a discontinuous motif once per CDR3
    combos <- unique(combos)
    return(combos)
  })
  return(discs)
}

# Get all of the discontinuous 5mers for the Cd154-
disc.fivemer.gen  <- function(df, CDR3.col){
  # Get rid of CDR3s that are too short to contribute
  CDR3.df <- df[nchar(df[[CDR3.col]]) >= 11,]
  
  # create all discontinuous nmers for the dataframe
  discs <- sapply(CDR3.df[[CDR3.col]], function(x){
    # Trim the CDR3
    trimmed.CDR3 <- substr(x, 4, nchar(x) - 3)
    
    # Make an empty vector
    combos <- c()
    # Loop across a CDR3 and get all possible discontinuous 3mers
    for(i in 1:(nchar(trimmed.CDR3) - 4)){
      comb.1 <- paste(substr(trimmed.CDR3, i, i), "X", substr(trimmed.CDR3, i + 2, i + 4), sep = "")
      comb.2 <- paste(substr(trimmed.CDR3, i, i + 1), "X", substr(trimmed.CDR3, i + 3, i + 4), sep = "")
      comb.3 <- paste(substr(trimmed.CDR3, i, i + 2), "X", substr(trimmed.CDR3, i + 4, i + 4), sep = "")
      comb.4 <- paste(substr(trimmed.CDR3, i, i), "XX", substr(trimmed.CDR3, i + 3, i + 4), sep = "")
      comb.5 <- paste(substr(trimmed.CDR3, i, i + 1), "XX", substr(trimmed.CDR3, i + 4, i + 4), sep = "")
      comb.6 <- paste(substr(trimmed.CDR3, i, i), "X", substr(trimmed.CDR3, i + 2, i + 2), "X",
                      substr(trimmed.CDR3, i + 4, i + 4), sep = "")
      combos <- c(combos, comb.1, comb.2, comb.3, comb.4, comb.5, comb.6)
    }
    # Only count a discontinuous motif once per CDR3
    combos <- unique(combos)
    return(combos)
  })
  return(discs)
}


# Function that will return CDR3s with a particular discontinuous nmer
find_disc_4mer <- function(nmer,df, CDR3.col, count.col, disc.list){
  ### Want to call disc.fourmer.gen here###
  
  # Create an empty dataframe that will be returned
  df.out <- data.frame(nmer = nmer,
                       unique.CDR3s = numeric(length = length(nmer)),
                       count = numeric(length = length(nmer)),
                       prop = numeric(length = length(nmer)))
  
  # Get rid of CDR3s that are too short to contribute
  CDR3.df <- df[nchar(df[[CDR3.col]]) >= 10,]
  
  list <- lapply(nmer, function(nmer){
    CDR3.list <- c()
    for(i in 1:length(disc.list)){
      if (nmer %in% disc.list[[i]]){
        CDR3.list <- c(CDR3.list, names(disc.list[i]))
      }
    }
    motif.df <- CDR3.df[CDR3.df[[CDR3.col]] %in% CDR3.list,]
    info <- list(length(unique(motif.df[[CDR3.col]])), sum(motif.df[[count.col]]),
                 sum(motif.df[[count.col]]) / sum(df[[count.col]]))
    return(info)
  })
  
  df.out$unique.CDR3s <- sapply(list, '[[', 1)
  df.out$count <- sapply(list, '[[', 2)
  df.out$prop <- sapply(list, '[[', 3)
  df.out;
  return(df.out)
}

# Function that will return CDR3s with a particular discontinuous nmer
find_disc_5mer <- function(nmer,df, CDR3.col, count.col, disc.list){
  ### Want to call disc.fourmer.gen here###
  
  # Create an empty dataframe that will be returned
  df.out <- data.frame(nmer = nmer,
                       unique.CDR3s = numeric(length = length(nmer)),
                       count = numeric(length = length(nmer)),
                       prop = numeric(length = length(nmer)))
  
  # Get rid of CDR3s that are too short to contribute
  CDR3.df <- df[nchar(df[[CDR3.col]]) >= 11,]
  
  list <- lapply(nmer, function(nmer){
    CDR3.list <- c()
    for(i in 1:length(disc.list)){
      if (nmer %in% disc.list[[i]]){
        CDR3.list <- c(CDR3.list, names(disc.list[i]))
      }
    }
    motif.df <- CDR3.df[CDR3.df[[CDR3.col]] %in% CDR3.list,]
    info <- list(length(unique(motif.df[[CDR3.col]])), sum(motif.df[[count.col]]),
                 sum(motif.df[[count.col]]) / sum(df[[count.col]]))
    return(info)
  })
  
  df.out$unique.CDR3s <- sapply(list, '[[', 1)
  df.out$count <- sapply(list, '[[', 2)
  df.out$prop <- sapply(list, '[[', 3)
  df.out;
  return(df.out)
}
