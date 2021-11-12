#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

base_dir        <- "/home/amy/GitHub/synthetic-data-genomics"
r_base          <- "research_paper_code"
experiment_dir  <- "mice_data_set" 
real_gwas_path = paste(base_dir, "/mice_data_set/out", sep="")
setwd(base_dir)

map <- function(){
    
    # Map QTLs for phenotypes measured in CFW outbred mice using the linear
    # mixed model (LMM) analysis implemented in GEMMA.

    # Set up working directory structures

    library(qtl)
    library(data.table)
    library(qqman)
    source(paste(r_base, "/src/misc.R", sep=""))
    source(paste(r_base, "/src/gemma.R", sep=""))
    source(paste(r_base, "/src/read.data.R", sep=""))
    source(paste(r_base, "/src/data.manip.R", sep=""))
    source(paste(r_base, "/src/qtl.analyses.R", sep=""))
    analysis_selection = analyses["abBMD"]
    for (analysis in analysis_selection) {
        print(analysis)
    }


    # SCRIPT PARAMETERS
    # -----------------
    chromosomes    <- NULL
    gemmadir       <- paste(experiment_dir, "/gemma", sep="")
    gemma.exe      <- paste("./", "gemma-0.98.4-linux-static-AMD64", sep="")
    geno_txt_base       <- paste(experiment_dir, "/data/synthetic_genomes", sep="")
    map_txt_base       <- paste(experiment_dir, "/data/map_abBMD", sep="")
    
    # Read in the synthetic phenotype data
    pheno_synth_file <- paste(experiment_dir, "/data/phenome_alldata_synth.csv", sep="")
    pheno_all <- read.csv(pheno_synth_file,quote = "",header = TRUE,check.names = FALSE,
                    stringsAsFactors = FALSE,comment.char = "#")
    
    geno_txt       <- paste(geno_txt_base,  ".txt", sep="")
    map_txt       <- paste(map_txt_base, ".txt", sep="")
    dt = run_gwas(geno_txt, map_txt, pheno_all, analysis, gemmadir, gemma.exe)
    dt
}    


run_gwas <- function(geno_txt, map_txt, pheno_all, analysis, gemmadir, gemma.exe) {

    # LOAD GENOTYPE DATA
    # ------------------
    # Load the "mean genotypes"; i.e., the the mean alternative allele
    # counts.

    map     <- read.map(map_txt)
    out     <- read.geno.dosage(geno_txt,nrow(map))
    discard <- out$discard
    X_all   <- out$geno
    rm(out)

    # Discard genotype samples from mislabeled flowcell samples.
    X_all <- X_all[which(discard == "no"),]


    # Discard SNPs with low "imputation quality" assessed by inspecting
    # the genotype probabilities. Retain SNPs for which: (1) at least 95%
    # of the samples have a maximum probability genotype greater than than
    # 0.5; (2) the minor allele frequency is greater than 2%.
    f       <- apply(X_all,2,compute.maf)
    markers <- which(map$quality > 0.95 & f > 0.02)
    map     <- map[markers,]
    X_all   <- X_all[,markers]
    
    # min_var, max_var - which of the columns in genotype data (X_all above) to be considered when 
    # using linear models (mostly for speed). Please note that gemma analysis will analyze the 
    # whole chromosome
    min_var = 1
    max_var = dim(X_all)[2]

    chromosomes <- unique(map[min_var:max_var,"chr"])
    
    ##################################
    # Cleanup data
    phenotype  <- analysis$pheno
    covariates <- analysis$cov
    outliers   <- analysis$outliers
 
    # Not sure if need
    pheno <- copy(pheno_all)
    if (!is.null(outliers))
      pheno_all <- remove.outliers(pheno,phenotype,covariates,outliers)

    
    # Only analyze samples (i.e. rows of the genotype and phenotype
    # matrices) for which the phenotype and all the covariates are
    # observed.
    pheno <- pheno[which(none.missing.row(pheno[c(phenotype,covariates)])),]
    
    # End not sure if need  
     
    # Align the phenotypes and genotypes
    ids   <- intersect(pheno_all$id,rownames(X_all))
    pheno <- pheno_all[match(ids,pheno_all$id),]
    X     <- X_all[match(ids,rownames(X_all)),]

    #Compute linear model

    #lm_out_csv <- paste(experiment_dir, "/out_synth/lm_", analysis$pheno, "_", min_var, "_", max_var, ".csv",sep="")
    lm_out_csv <- paste(experiment_dir, "/out_synth/lm_", analysis$pheno, ".csv",sep="")
    
    if(!file.exists(lm_out_csv)) {
        print(dim(X)[2])
        dt <- data.table(snp=rep("",dim(X)[2]), chr=rep(0,dim(X)[2]), pos=rep(0,dim(X)[2]), p=rep(1,dim(X)[2]))
        for (i in min_var:max_var) {
            X_variant <- cbind(X[,i], pheno_column=pheno[,analysis$pheno])
            colnames(X_variant)[1]<-colnames(X)[i]
            f <- paste("pheno_column ~ ",colnames(X)[i])
            # Add any covariates
            for(cov in analysis$cov) {
                X_variant <- cbind(X_variant, pheno_column=pheno[,cov])
                f <- paste(f,"+",cov)
            }
            res_variant <- lm(pheno_column~., data = data.frame(X_variant))
          
            dt[i,1] = colnames(X)[i]
            dt[i,2] = as.numeric(map[map["id"]==colnames(X)[i],"chr"])
            dt[i,3] = as.numeric(map[map["id"]==colnames(X)[i],"pos"])
            dt[i,4] = as.numeric(summary(res_variant)$coefficients[2,4])
        }

        # Print to file
        #write.csv(dt[order(rank(p)),][1:(max_var-min_var)], lm_out_csv)
        write.csv(dt, lm_out_csv)

    }
    
    dt
    
}

dt = map()
