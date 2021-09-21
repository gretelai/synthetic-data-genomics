# SUMMARY
# -------
# This file contains some function definitions that don't fit in any
# other file. Here is an overview of the functions defined in this
# file:
#
#    caterase(s)
#    logit10(x)
#    project.onto.interval(x,a,b)
#    center.columns(X)
#    is.binary.factor(x)
#    is.binary(x)
#    binfactor2num(x)
#    binary.from.categorical(x,col.names)
#    none.missing.row(x)
#    compute.maf(geno)
#    check.normal.quantiles(x)
#
# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------
# Output the string using 'cat', then move the cursor back to the
# beginning of the string so that subsequent output will overwrite
# this string.
caterase <- function (s)
  cat(s,rep("\b",nchar(s)),sep="")

# ----------------------------------------------------------------------
# Return the logit of x. (Here we use the base-10 logarithm instead of the 
# more commonly used natural logarithm.)
logit10 <- function (x) {
  eps <- .Machine$double.eps
  return(log10((x + eps)/(1 - x + eps)))
}

# ----------------------------------------------------------------------
# Return x projected onto interval [a,b].
project.onto.interval <- function (x, a, b)
  pmin(b,pmax(a,x))

# ----------------------------------------------------------------------
# Center the columns of matrix X so that the entries in each column
# of X add up to 0.
# 
center.columns <- function (X) {
  y <- colMeans(X)
  X <- X - matrix(y,nrow(X),ncol(X),byrow = TRUE)
  return(X)
}

# ----------------------------------------------------------------------
# Return TRUE if x is a factor with exactly 2 levels.
is.binary.factor <- function (x)
  is.factor(x) & nlevels(x) == 2
    
# ----------------------------------------------------------------------
# Return TRUE if all elements of a numeric vector are either 0 or 1,
# ignoring missing values (NA).
is.binary <- function (x)
  is.numeric(x) & sum(!(x == 0 | x == 1),na.rm = TRUE) == 0

# ----------------------------------------------------------------------
# Convert a factor with exactly 2 levels to a numeric vector with
# values 0 or 1.
binfactor2num <- function (x) {
  if (!is.binary.factor(x))
    stop("Factor must have exactly 2 levels")
  return(as.numeric(x) - 1)
}
  
# ----------------------------------------------------------------------
# Return a data frame with one column for each level of factor x. The
# columns of the data frame encode the categorical variable with n
# levels as n binary variables.
binary.from.categorical <- function (x, col.names = NULL) {

  # Create a binary factor for each value of the categorical variable.
  d <- list()
  for (value in levels(x))
    d[[value]] <- factor(as.integer(x == value))

  # Convert the list to a data frame, and adjust the column names, if
  # requested.
  d <- data.frame(d,check.names = FALSE)
  if (!is.null(col.names))
    names(d) <- col.names

  # Output the newly created data frame.
  return(d)
}

# ----------------------------------------------------------------------
# For each row of the matrix or data frame, return TRUE if all the
# entries in the row are provided (that is, they are not missing).
none.missing.row <- function (x)
  rowSums(is.na(x)) == 0

# ----------------------------------------------------------------------
# Returns the minor allele frequency given a vector of genotypes
# encoded as allele counts.
compute.maf <- function (geno) {
  f <- mean(geno,na.rm = TRUE)/2
  return(min(f,1-f))
}

# ----------------------------------------------------------------------
# Compare the empirical quantiles of x against the expected valuesx
# under the normal distribution.
check.normal.quantiles <- function (x) {

  # Discard the missing values.
  x <- x[!is.na(x)]
    
  # Transform the observations so that they have zero mean and unit
  # variance.
  x <- (x - mean(x))/sqrt(var(x))
  
  # Create a data frame giving the observed and expected proportion of
  # samples within 1, 2 and 3 standard deviations of the mean.
  return(data.frame(exp = c(pnorm(1) - pnorm(-1),
                            pnorm(2) - pnorm(-2),
                            pnorm(3) - pnorm(-3)),
                    obs = c(mean(-1 < x & x < 1),
                            mean(-2 < x & x < 2),
                            mean(-3 < x & x < 3)),
                    row.names = c("sd1","sd2","sd3")))
}
