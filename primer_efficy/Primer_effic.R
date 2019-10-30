#!/usr/bin/env Rscript
# Script to caculate primer efficeincy 
# Assumes a df of at least three columns:
#                                    Sample.Name : cDNA dilution "1:x"
#                                    CT          : qPCR CT value
#                                    Target.Name : Name of gene "primer_pair"
# Useage: Rscript Primer_effic.R ~/Desktop/prime_eff.csv h2o :
# Rscript Primer_effic.R ~/Desktop/B1_PE.csv h2o
# 
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) == 2)
  x <- read.csv(args[1])
  # These names correspond to the output of a specific machine
  dat <- x[c("Sample.Name", "CT", "Target.Name")]
  # Set the id of the Target name to simplify things in the feature
  # Process
  primer <- unique(dat["Target.Name"])
  #
  forma <- sapply(primer,
                  function(x) paste(as.character(x)))
  
  
  # Target.Name is the primer name
  umods <- lapply(forma, function(x){subset(dat, Target.Name == x)})
  
  # Prepare aggregation function
  
  # Function to force numeric
  mean_convert <- function(x) {
    mean(as.numeric(as.character(x)))
  }
  
  # Function to subset only numeric variables
  numeric_columns <- function(x){
    num_col <- sapply(x, is.numeric)
    return(x[ , num_col])
  }
  
  # Function to aggregate CT values by mean
  # You may Use 
  aggregate_CT <- function(x) {
    # Configure for quntstudio's output "Sample.Name, Target.Name"
    x_num <- numeric_columns(x)
    aggregated <- aggregate(x_num, by = list(x$Sample.Name),
                            FUN = mean_convert)
    names(aggregated) <- c("Sample_", "CT_mean")
    aggregated$Target <- unique(x$Target.Name)
    return(aggregated)
  }
  
  # Aggregated data
  agg_data <- lapply(umods, aggregate_CT)
  
  # Function to format the dultion
  create_dilute <- function(x){
    xx <- subset(x, !Sample_ == args[2])
    xx$s_Dil <- log((as.numeric(xx$Sample_)))
    return(xx)
  }
  
  #
  ad_dilute <- lapply(agg_data, create_dilute)
  
  # Calculate efficiency
  res <- lapply(ad_dilute,
                function(x){s <- lm(CT_mean ~ s_Dil, x)$coefficients[2]
                s <- 10^(-1/s)
                s <- (s - 1)*100
                s <- paste(unique(x$Target), s, sep = " : ")
                })
  # Print results
  for (i in res) {
    cat(i, sep = "\n")
  }
}

main()
