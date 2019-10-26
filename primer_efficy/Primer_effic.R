# Function to caculate primer efficeincy 
# Assumes a df of least three columns:
#                                    Sample.Name : cDNA dilution "1:x"
#                                    CT          : qPCR CT value
#                                    Target.Name : Name of gene "primer_pair"
primer_efficy <- function(x,H2O_name, sep_dilution) {
  # These names correspond to out of a specific machine
  # Check your output
  #
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
  
  # Function to format the dultion to plot able
  create_dilute <- function(x){
    xx <- subset(x, !Sample_ == H2O_name)
    xx$Dil <- gsub(paste("*.", sep_dilution, sep = ""), "", xx$Sample_)
    xx$s_Dil <- log(1/as.numeric(xx$Dil))
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
    print(i)
  }
}

# Calculate example
foo <- read.csv("~/Desktop/prime_eff.csv")
primer_efficy(foo, "h2o", ":")
