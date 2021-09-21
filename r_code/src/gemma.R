# SUMMARY
# -------
# This file defines some functions that I use for QTL mapping in
# GEMMA. Here is an overview of the functions defined in this file:
#
#   write.gemma.pheno(file,phenotype,pheno)
#   write.gemma.covariates(file,covariates,pheno) 
#   write.gemma.map(file,map)
#   write.gemma.geno(file,geno,map)
#   read.gemma.assoc(file)
#   run.gemma(phenotype,covariates,pheno,geno,map,gemmadir,gemma.exe,
#             chromosomes)
#   run.gemma.norr(phenotype,covariates,pheno,geno,map,gemmadir,gemma.exe) 
#
# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------
# Write the phenotype data to a file in the format used by GEMMA. Each
# line of the file contains one phenotype observation.
write.gemma.pheno <- function (file, phenotype, pheno) {
  y <- pheno[phenotype]
  if (is.numeric(y))
    y <- round(y,digits = 6)
  write.table(y,file,quote = FALSE,row.names = FALSE,col.names = FALSE)
}

# ----------------------------------------------------------------------
# Write the covariate data to a file in the format used by GEMMA. Each
# line corresponds to a sample. We must include an additional
# covariate for the intercept.
write.gemma.covariates <- function (file, covariates, pheno) {
  if (is.null(covariates)) {
    write.table(data.frame(rep(1,nrow(pheno))),file,sep = " ",
                quote = FALSE,row.names = FALSE,col.names = FALSE)
  } else {
    write.table(cbind(1,data.frame(lapply(pheno[covariates],function (x) {
      if (is.numeric(x))
        round(x,digits=6)
      else
        x
    }))),file,sep = " ",quote = FALSE,row.names = FALSE,col.names = FALSE)
  }
}

# ----------------------------------------------------------------------
# Write the SNP information to a space-delimited text file in the
# format used by GEMMA. This file contains one line per SNP, with
# three columns: (1) SNP label, (2) base-pair position, (3)
# chromosome.
write.gemma.map <- function (file, map)
  write.table(map[c("id","pos","chr")],file,sep = " ",quote = FALSE,
              row.names = FALSE,col.names = FALSE)

# ----------------------------------------------------------------------
# Store the mean genotypes as a space-delimited text file in the
# format used by GEMMA, in which we have one row per SNP, and one
# column per sample. The first three column give the SNP label, and
# the two alleles.
write.gemma.geno <- function (file, geno, map) {
  geno <- t(geno)
  geno <- as.data.frame(geno,check.names = FALSE)
  geno <- round(geno,digits = 3)
  geno <- cbind(map[c("id","ref","alt")],geno)
  write.table(geno,file,sep = " ",quote = FALSE,row.names = FALSE,
              col.names = FALSE)
}

# ----------------------------------------------------------------------
# Reads in the association results from GEMMA, and returns a data
# frame containing three columns: chromosome number ("chr"); base-pair
# position ("pos"); and the base-10 logarithm of the p-value ("log10p").
read.gemma.assoc <- function (file) {
  gwscan <- read.table(file,sep = "\t",header = TRUE,check.names = FALSE,
                  quote = "",stringsAsFactors = FALSE)
  rownames(gwscan) <- gwscan$rs
  gwscan           <- gwscan[c("chr","ps","p_lrt")]
  gwscan           <- transform(gwscan,p_lrt = -log10(p_lrt))
  colnames(gwscan) <- c("chr","pos","log10p")
  return(gwscan)
}

# ----------------------------------------------------------------------
# This function maps QTLs using GEMMA, writing all the files required
# by GEMMA to the directory specified by "gemmadir". The QTLs are
# mapped separately for each chromosome, in which the kinship matrix
# is computed using all markers except the markers on the given
# chromosome.
run.gemma <- function (phenotype, covariates, pheno, geno, map,
                       gemmadir, gemma.exe, chromosomes = NULL) {

  # Give summary of analysis.
  cat("Mapping QTLs for",phenotype,"in",nrow(pheno),"mice, ")
  if (!is.null(covariates)) {
    cat("controlling for ",paste(covariates,collapse=" + "),".\n",sep="")
  } else {
    cat("with no covariates included.\n")
  }

  # Write the phenotype and covariate data to separate files.
  cat("Writing phenotype and covariate data to file.\n")
  write.gemma.pheno(paste0(gemmadir,"/pheno.txt"),phenotype,pheno)
  write.gemma.covariates(paste0(gemmadir,"/covariates.txt"),covariates,pheno) 

  # We will map QTLs on these chromosomes.
  if (is.null(chromosomes))  
    chromosomes  <- levels(map$chr)
  
  # Repeat for each chromosome.
  scans        <- vector("list",length(chromosomes))
  names(scans) <- chromosomes

  for (chr in chromosomes) {
    # Compute the kinship matrix using all markers that are *not* on
    # the chromosome.
    cat("Mapping QTLs on chromosome ",chr,".\n",sep="")
    cat(" * Computing kinship matrix.\n")
    markers <- which(map$chr != chr)
    K <- tcrossprod(center.columns(geno[,markers])) / length(markers)

    # Save the kinship matrix to a text file.
    cat(" * Writing kinship matrix to file.\n")
    write.table(round(K,digits = 6),paste0(gemmadir,"/kinship.txt"),sep = " ",
                quote = FALSE,row.names = FALSE,col.names = FALSE)
  
    # Write out the mean genotypes and map information for all markers
    # on the chromosome.
    markers <- which(map$chr == chr)
    cat(" * Writing to file genotypes for ",length(markers),
        " markers on chromosome ",chr,".\n",sep="")
    write.gemma.geno(paste0(gemmadir,"/geno.txt"),geno[,markers],map[markers,])
    cat(" * Writing genetic map for",length(markers),
        "markers on chromosome",chr,"to file.\n")
    write.gemma.map(paste0(gemmadir,"/map.txt"),map[markers,])    

      
    # Set the local directory to the location of the GEMMA files.
    srcdir <- getwd()
    setwd(gemmadir)
      
    # Now we are finally ready to run GEMMA for all markers on the
    # chromosome using the kinship matrix computed using all the
    # markers *not* on the chromosome.
    cat(" * Computing p-values for ",length(markers),
        " markers on chromosome ",chr,".\n",sep="")
    system(paste(gemma.exe,"-g geno.txt -a map.txt -p pheno.txt",
                 "-c covariates.txt -k kinship.txt -notsnp -lmm 2",
                 "-lmin 0.01 -lmax 100"),
           ignore.stdout = TRUE)
      
    # Restore the working directory.
    setwd(srcdir)
      
    # Load the results of the GEMMA association analysis.
    scans[[chr]] <-
      read.gemma.assoc(paste0(gemmadir,"/output/result.assoc.txt"))
  }


  
  # Merge the mapping results from all chromosomes into a single table.
  gwscan           <- do.call(rbind,scans)
  rownames(gwscan) <- do.call(c,lapply(scans,rownames))
  class(gwscan)    <- c("scanone","data.frame")

  # Return the genome-wide scan.
  return(gwscan)
}

# ----------------------------------------------------------------------
# This function maps QTLs using GEMMA without the "realized
# relatedness" matrix to account for population structure.
run.gemma.norr <- function (phenotype, covariates, pheno, geno, map,
                            gemmadir, gemma.exe) {


  # Give summary of analysis.
  cat("Mapping QTLs for",phenotype,"in",nrow(pheno),"mice, ")
  if (!is.null(covariates)) {
    cat("controlling for ",paste(covariates,collapse=" + "),".\n",sep="")
  } else {
    cat("with no covariates included.\n")
  }

  # Write the phenotype and covariate data to separate files.
  cat("Writing phenotype and covariate data to file.\n")
  write.gemma.pheno(paste0(gemmadir,"/pheno.txt"),phenotype,pheno)
  write.gemma.covariates(paste0(gemmadir,"/covariates.txt"),covariates,pheno) 

  # Write out the mean genotypes and map information for all markers.
  cat("Writing SNP and genotype data to file.\n")
  write.gemma.map(paste0(gemmadir,"/map.txt"),map)
  write.gemma.geno(paste0(gemmadir,"/geno.txt"),geno,map)
  
  # Write out the kinship matrix to file.
  cat("Writing identity kinship matrix to file.\n");
  write.table(diag(nrow(pheno)),paste0(gemmadir,"/kinship.txt"),
              sep = " ",quote = FALSE,row.names = FALSE,
              col.names = FALSE)
  # cat("Writing full kinship matrix to file.\n");
  # K <- tcrossprod(center.columns(geno)) / nrow(map)
  # write.table(round(K,digits = 6),paste0(gemmadir,"/kinship.txt"),
  #             sep = " ",quote = FALSE,row.names = FALSE,
  #             col.names = FALSE)

  # Set the local directory to the location of the GEMMA files.
  srcdir <- getwd()
  setwd(gemmadir)

  # Now we are finally ready to run GEMMA for all markers.
  cat("Computing p-values for",nrow(map),"candidate markers.\n")
  system(paste(gemma.exe,"-g geno.txt -a map.txt -p pheno.txt",
               "-c covariates.txt -k kinship.txt -notsnp -lmm 2",
               "-lmin 0.01 -lmax 100"),
         ignore.stdout = TRUE)

  # Restore the working directory.
  setwd(srcdir)

  # Load the results of the GEMMA association analysis.
  gwscan <- read.gemma.assoc(paste0(gemmadir,"/output/result.assoc.txt"))
  class(gwscan) <- c("scanone","data.frame")

  # Restore the working directory.
  setwd(srcdir)

  # Return the genome-wide scan.
  return(gwscan)
}
