# SUMMARY
# -------
# This file contains several functions for reading the QTL experiment
# data from text files. Here is an overview of the functions defined
# in this file:
#
#   read.pheno(file)
#   read.map(file)
#   read.geno.dosage(file,n)
#
# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------
# Loads the phenotype data stored in a CSV file.
read.pheno <- function (file) {
    
  # Read in the phenotype data from the CSV file. 
  pheno <- read.csv(file,quote = "",header = TRUE,check.names = FALSE,
                    stringsAsFactors = FALSE,comment.char = "#")

  # Convert some of the columns to factors.
  pheno <- transform(pheno,
                     id            = as.character(id),
                     round         = factor(round,paste0("SW",1:25)),
                     FCbox         = factor(FCbox),
                     PPIbox        = factor(PPIbox),
                     methcage      = factor(methcage),
                     methcycle     = factor(methcycle),
                     discard       = factor(discard),
                     mixup         = factor(mixup),
                     earpunch      = factor(earpunch),
                     abnormalbone  = factor(abnormalbone),
                     experimenters = factor(experimenters))

  # Convert some of the columns to double precision.
  pheno <- transform(pheno,
                     fastglucose       = as.double(fastglucose),
                     D1totaldist0to15  = as.double(D1totaldist0to15),
                     D1totaldist15to30 = as.double(D1totaldist15to30),
                     D1totaldist0to30  = as.double(D1totaldist0to30),
                     D2totaldist0to15  = as.double(D2totaldist0to15),
                     D2totaldist15to30 = as.double(D2totaldist15to30),
                     D2totaldist0to30  = as.double(D2totaldist0to30),
                     D3totaldist0to15  = as.double(D3totaldist0to15),
                     D3totaldist15to30 = as.double(D3totaldist15to30),
                     D3totaldist0to30  = as.double(D3totaldist0to30),

                     D1TOTDIST5  = as.double(D1TOTDIST5),
                     D1TOTDIST10 = as.double(D1TOTDIST10),
                     D1TOTDIST15 = as.double(D1TOTDIST15),
                     D1TOTDIST20 = as.double(D1TOTDIST20),
                     D1TOTDIST25 = as.double(D1TOTDIST25),
                     D1TOTDIST30 = as.double(D1TOTDIST30),

                     D2TOTDIST5  = as.double(D2TOTDIST5),
                     D2TOTDIST10 = as.double(D2TOTDIST10),
                     D2TOTDIST15 = as.double(D2TOTDIST15),
                     D2TOTDIST20 = as.double(D2TOTDIST20),
                     D2TOTDIST25 = as.double(D2TOTDIST25),
                     D2TOTDIST30 = as.double(D2TOTDIST30),

                     D3TOTDIST5  = as.double(D3TOTDIST5),
                     D3TOTDIST10 = as.double(D3TOTDIST10),
                     D3TOTDIST15 = as.double(D3TOTDIST15),
                     D3TOTDIST20 = as.double(D3TOTDIST20),
                     D3TOTDIST25 = as.double(D3TOTDIST25),
                     D3TOTDIST30 = as.double(D3TOTDIST30))

  # Return the phenotype data table.
  return(pheno)
}

# ----------------------------------------------------------------------
# Returns a data frame containing the marker data stored in a text
# file. Here I convert the chromosomes and alleles to factors manually
# to make sure that the chromosomes and bases are ordered properly in
# the factors.
read.map <- function (file) {
  bases <- c("A","T","G","C")
  map   <- read.table(file,sep = " ",header = TRUE,comment.char = "#",
                 stringsAsFactors = FALSE)
  return(transform(map,
                   chr = factor(chr,1:19),
                   ref = factor(ref,bases),
                   alt = factor(alt,bases)))
}

# ----------------------------------------------------------------------
# Returns a list object containing (1) a vector "discard" specifying
# samples that are from potentially mislabeled flow samples, (2) an n
# x p matrix of genotype "dosages" (expected allele counts), where n
# is the number of samples and p is the number of SNPs.
read.geno.dosage <- function (file, n) {
  out <- fread(file,sep = " ",header = TRUE,showProgress = FALSE,
               colClasses = c("character","character",rep("double",n)))
  class(out)    <- "data.frame"
  rownames(out) <- out$id
  return(list(discard = factor(out$discard),
              geno    = as.matrix(out[-(1:2)])))
}
