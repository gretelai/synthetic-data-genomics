# SUMMARY
# -------
# Some functions for manipulating QTL experiment data. Here is an
# overview of the functions defined in this file:
#
#    prepare.pheno(pheno)
#    remove.outliers(pheno,phenotype,covariates,outliers,verbose)
#    genotypes2counts(geno,A,B)
#    mapgeno(gp)
#    data2qtl(pheno,gp,map)
#    merge.gwscan.qtl(files)
#    get.thresholds(perms,alpha)
#
# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------
# This function prepares the phenotype data for QTL mapping.
prepare.pheno <- function (pheno) {

  # Remove rows marked as "discard" and as possible sample mixups.
  pheno <- subset(pheno,discard == "no" & mixup == "no")

  # Create a new binary trait that is equal to 1 when the mouse has an
  # "abnormal" bone-mineral density, and 0 otherwise.
  pheno <- cbind(pheno,
                 data.frame(abBMD = factor(as.numeric(pheno$BMD > 90))))

  # Create some binary traits ("pretrainbin", "d1tonebin" and
  # "altcontextbin") for fear conditioning metrics that show a
  # "bimodal" freezing distribution. Each of these phenotypes is equal
  # to 1 when we observe some freezing, and 0 otherwise.
  pheno <-
    cbind(pheno,with(pheno,
        data.frame(pretrainbin   = factor(as.numeric(PreTrainD1     > 0.01)),
                   d1tonebin     = factor(as.numeric(AvToneD1       > 0.01)),
                   altcontextbin = factor(as.numeric(AvAltContextD3 > 0.01)))))

  # Create a binary trait from the p120b1 trait that indicates that
  # the mouse does not respond to the pulse, possibly indicating
  # deafness. I choose mice with p120b1 values on the "long tail" of
  # the distribution for this phenotype.
  pheno <- cbind(pheno,
                 data.frame(deaf = factor(as.numeric(pheno$p120b1 < 50))))

  # The deaf mice effectively add "noise" to the PPI phenotypes, in
  # part because the ratio has a denominator that would be close to
  # zero. So I remove these samples from the analysis.
  PPI.phenotypes <-
      c("p120b1","p120b2","p120b3","p120b4","pp3b1","pp3b2","pp3avg",
        "pp6b1","pp6b2","pp6avg","pp12b1","pp12b2","pp12avg","pp3PPIb1",
        "pp3PPIb2","pp3PPIavg","pp6PPIb1","pp6PPIb2","pp6PPIavg",
        "pp12PPIb1","pp12PPIb2","pp12PPIavg","PPIavg","startle")
  rows <- which(pheno$p120b1 < 10)
  pheno[rows,PPI.phenotypes] <- NA
  
  # Transform bone-mineral densities, which are ratios, to the
  # log-scale (in base 10), so that the distribution of this phenotype
  # looks roughly normal. (Although there is a long tail of high BMD
  # measurements that we will have to discard.)
  pheno <- transform(pheno,BMD = log10(BMD))

  # Transform the fear conditioning phenotypes, which are proportions
  # (numbers between 0 and 1) to the log-odds scale using the logit
  # function (in base 10).
  f     <- function (x) logit10(project.onto.interval(x,0.01,0.99))
  pheno <- transform(pheno,
                     PreTrainD1     = f(PreTrainD1),
                     AvToneD1       = f(AvToneD1),
                     AvContextD2    = f(AvContextD2),
                     AvAltContextD3 = f(AvAltContextD3),
                     AvToneD3       = f(AvToneD3),
                     D3.180         = f(D3.180),
                     D3.240         = f(D3.240),
                     D3.300         = f(D3.300),
                     D3.360         = f(D3.360))

  # Transform some of the prepulse inhibition (PPI) phenotypes, which
  # are proportions (numbers between 0 and 1, except for a few
  # negative values), to the log-odds scale.
  pheno <- transform(pheno,
                     pp3PPIavg  = f(pp3PPIavg),
                     pp6PPIavg  = f(pp6PPIavg),
                     pp12PPIavg = f(pp12PPIavg),
                     PPIavg     = f(PPIavg))

  # Transform the "center time" phenotypes so that they are
  # proportions, then transform the proportions to the (base 10) logit
  # scale.
  pheno <- transform(pheno,
                     D1ctrtime0to15 = f(D1ctrtime0to15/900),
                     D2ctrtime0to15 = f(D2ctrtime0to15/900),
                     D3ctrtime0to15 = f(D3ctrtime0to15/900),
                     
                     D1ctrtime0to30 = f(D1ctrtime0to30/1800),
                     D2ctrtime0to30 = f(D2ctrtime0to30/1800),
                     D3ctrtime0to30 = f(D3ctrtime0to30/1800))

  # Transform the "vertical activity" phenotypes so that they are
  # proportions, then transform the proportions to the (base 10) logit
  # scale.
  pheno <- transform(pheno,
                     D1vact0to15 = f(D1vact0to15/900),
                     D2vact0to15 = f(D2vact0to15/900),
                     D3vact0to15 = f(D3vact0to15/900),

                     D1vact0to30 = f(D1vact0to30/1800),
                     D2vact0to30 = f(D2vact0to30/1800),
                     D3vact0to30 = f(D3vact0to30/1800))
                     
  # Create a new phenotype defined to be the first principal component
  # of the muscle weights.
  mw.pheno           <- pheno[c("TA","EDL","soleus","plantaris","gastroc")]
  rows               <- which(rowSums(is.na(mw.pheno)) == 0)
  pheno              <- cbind(pheno,data.frame(mwpc = NA))
  pheno[rows,"mwpc"] <- prcomp(mw.pheno[rows,],center=TRUE,scale=TRUE)$x[,1]

  # Return the updated phenotype data table.
  return(pheno)
}

# ----------------------------------------------------------------------
# Remove outliers from the phenotype data for the given phenotype,
# optionally conditioning on covariates. If covariates are specified,
# outliers are determined according to the residual of the phenotype
# after regressing on the covariates.
remove.outliers <- function (pheno, phenotype, covariates, outliers,
                             verbose = TRUE) {

  # Specify the function for removing the outliers.
  is.outlier <- function (x) {
    y           <- outliers(x)
    y[is.na(y)] <- FALSE
    return(y)
  }
  
  # If we are conditioning on one or more covariates, get the
  # residuals of the phenotype conditioned on these covariates.
  if (length(covariates) > 0) {
    f <- formula(paste(phenotype,"~",paste(covariates,collapse="+")))
    r <- resid(lm(f,pheno,na.action = na.exclude))
  } else
    r <- pheno[[phenotype]]

  # If requested, report how many outlying data points are removed.
  if (verbose) {
    n <- sum(is.outlier(r))
    if (n == 0)
      cat("No outliers for",phenotype)
    else
      cat(n,"outliers for",phenotype)
    if (length(covariates) == 0)
      cat(" are removed.\n")
    else
      cat(" conditioned on",paste(covariates,collapse=" + "),"are removed.\n")
  }
    
  # Remove the outlying data points.
  pheno[is.outlier(r),phenotype] <- NA

  # Return the updated phenotype data table.
  return(pheno)
}

# ----------------------------------------------------------------------
# This function takes three inputs: a vector of strings representing
# genotypes, a character giving the reference allele, and another
# character giving the alternative allele. The return value is a
# vector of allele counts; precisely, the number of times the
# alternative (B) allele appears in the genotype.
genotypes2counts <- function (geno, A, B) {

  # Initialize the allele counts to be missing.
  g <- rep(NA,length(geno))
    
  # Set the allele counts according to the genotypes.
  g[geno == paste0(A,A)] <- 0
  g[geno == paste0(A,B)] <- 1
  g[geno == paste0(B,A)] <- 1
  g[geno == paste0(B,B)] <- 2

  # Return the allele counts.
  return(g)
}

# ----------------------------------------------------------------------
# Given the genotype probabilities, output the genotypes (allele
# counts) with the highest probabilities; i.e. the maximum a
# posteriori estimates of the genotypes.
mapgeno <- function (gp) {

  # Get the largest probability for each genotype.
  r <- pmax(gp$p0,gp$p1,gp$p2)

  # Initialize the output.
  geno   <- r
  geno[] <- 0
  
  # Get the genotypes corresponding to the largest probabilities.
  geno[gp$p1 == r] <- 1
  geno[gp$p2 == r] <- 2

  # Return the matrix of genotypes.
  return(geno)
}

# ----------------------------------------------------------------------
# Convert the QTL experiment data to the format used in R/qtl. The
# return value is a 'cross' object that keeps track of all the data in
# a single QTL experiment; for more details, see the help for function
# 'read.cross' in R/qtl. The alleles are labeled as A and B, where A
# is the reference allele, and B is the alternative allele. The qtl
# library requires a paternal grandmother ("pgm") phenotype, so a
# column in the table for this phenotype is included if it is missing,
# and the entries of this column are set to zero. IMPORTANT NOTE: this
# function does not work for markers on the X chromosome.
data2qtl <- function (pheno, gp, map) {

  # Get the number of markers (p) and the number of samples (n).
  p <- nrow(map)
  n <- nrow(pheno)

  # Get the set of chromosomes.
  chromosomes <- as.character(unique(map$chr))
  
  # Create a new matrix in which the genotypes are encoded as follows:
  # NA = missing, 1 = AA, 2 = AB, 3 = BB. Here, A refers to the
  # reference allele, and B refers to the alternative allele (this
  # matches the encoding for the F2 intercross; see the help for
  # 'read.cross' in the qtl library for more details). Any genotypes
  # for which we are somewhat uncertain (probability less than 0.95),
  # I set these genotypes to missing (NA). 
  cutoff <- 0.95
  geno   <- matrix(NA,nrow = n,ncol = p,dimnames = list(pheno$id,map$snp))
  geno[gp$p0 > cutoff] <- 1
  geno[gp$p1 > cutoff] <- 2
  geno[gp$p2 > cutoff] <- 3

  # Convert the genotype probabilities to an n x p x 3 array.
  prob <- array(NA,c(n,p,3),list(pheno$id,map$snp,c("AA","AB","BB")))
  prob[,,1] <- gp$p0
  prob[,,2] <- gp$p1
  prob[,,3] <- gp$p2
  rm(gp)
  
  # Add the paternal grandmother ("pgm") phenotype if it does not
  # already exist.
  if (!is.element("pgm",names(pheno)))
    pheno <- cbind(pheno,pgm = 0)
  
  # Initialize the cross object, setting the genotypes to an empty list.
  cross        <- list(geno = list(),pheno = pheno)
  class(cross) <- c("f2","cross")

  # Set the alleles.
  attributes(cross)$alleles <- c("A","B")

  # Split up the markers by chromosome.
  for (i in chromosomes) {

    # Get the markers on the current chromosome.
    markers <- which(map$chr == i)

    # Get the genotype and map data corresponding to these markers.
    # Note that I need to coerce the genotype data to be a matrix in
    # the exceptional case when there is only a single marker
    # genotyped on the chromosome.
    m <- map[markers,]
    d <- list(data = as.matrix(geno[,markers]),
              prob = prob[,markers,],
              map  = m$pos)    
    attributes(d$prob)$map <- d$map
    names(d$map)           <- m$snp
    class(d)               <- "A"
    
    # Store the genotypes for markers on the chromosome.
    cross$geno[[i]] <- d    
  }

  return(cross)
}

# ----------------------------------------------------------------------
# Merge the results of R/qtl genome-wide scans for multiple phenotypes into
# a single data table.
merge.gwscan.qtl <- function (files) {

  # First load the marker data, and initialize the data table
  # containing the QTL mapping results.
  load(files[[1]])
  class(gwscan.qtl)       <- "data.frame"
  gwscan.merged           <- gwscan.qtl[c("chr","pos")]
  rownames(gwscan.merged) <- map$snp

  # Initialize other data structures containing the merged mapping
  # results.
  phenotypes.merged <- list()
  covariates.merged <- list()
  perms.merged      <- list()

  # Repeat for each file containing QTL mapping results.
  for (analysis in names(files)) {

    # Load the QTL mapping results.
    load(files[[analysis]])

    # Get the phenotype and covariates.
    phenotypes.merged[analysis] <- list(phenotype)
    covariates.merged[analysis] <- list(covariates)
    
    # Get the QTL mapping results.
    class(gwscan.qtl)    <- "data.frame"
    names(gwscan.qtl)[3] <- analysis
    gwscan.merged        <- cbind(gwscan.merged,gwscan.qtl[analysis])

    # Get the permutations.
    class(perms.qtl)         <- "matrix"
    perms.merged[[analysis]] <- perms.qtl
  }

  # Convert the permutations to a matrix, in which matrix columns
  # correspond to phenotypes.
  perms.merged           <- do.call(cbind,perms.merged)
  colnames(perms.merged) <- names(files)
  
  # Output the merged data.
  phenotypes.merged    <- unlist(phenotypes.merged)
  class(gwscan.merged) <- c("scanone","data.frame")
  class(perms.merged)  <- c("scanoneperm","matrix")
  return(list(phenotypes = phenotypes.merged,
              covariates = covariates.merged,
              gwscan.qtl = gwscan.merged,
              perms.qtl  = perms.merged))
}

# ----------------------------------------------------------------------
# Get the genome-wide LOD thresholds estimated in R/qtl using a
# permutation-based test.
get.thresholds <- function (perms, alpha) {
  thresholds        <- summary(perms,alpha)
  cols              <- colnames(thresholds)
  thresholds        <- as.vector(thresholds)
  names(thresholds) <- cols
  return(thresholds)
}
