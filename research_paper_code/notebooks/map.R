# Script version of map.ipynb
# To run, use the following shell command:
# R < map.R --no-save 2> r_stdout.txt 1> r_stderr.txt

# Set up working directory structures
library(stringr)
base_dir        <- str_replace(getwd(), 'research_paper_code/notebooks', '')
r_base          <- "research_paper_code"
experiment_dir  <- "mice_data_set"

setwd(base_dir)
getwd()

# Map QTLs for phenotypes measured in CFW outbred mice using the linear
# mixed model (LMM) analysis implemented in GEMMA.
library(qtl)
library(data.table)
source(paste(r_base, "/src/misc.R", sep=""))
source(paste(r_base, "/src/gemma.R", sep=""))
source(paste(r_base, "/src/read.data.R", sep=""))
source(paste(r_base, "/src/data.manip.R", sep=""))
source(paste(r_base, "/src/qtl.analyses.R", sep=""))

# SCRIPT PARAMETERS
# -----------------
chromosomes    <- NULL
gemmadir       <- paste(experiment_dir, "/gemma", sep="")
gemma.exe      <- paste("./", "gemma-0.98.4-linux-static-AMD64", sep="")
geno_txt       <- paste(experiment_dir, "/data/geno.txt", sep="")
map_txt        <- paste(experiment_dir, "/data/map.txt", sep="")
pheno_csv      <- paste(experiment_dir, "/data/pheno.csv", sep="")

# LOAD PHENOTYPE DATA
# -------------------
# Load the phenotype data, and discard outlying phenotype values. I
# create binary covariates from some of the categorical phenotypes.
cat("Loading phenotype data.\n")
pheno_all <- read.pheno(pheno_csv)
pheno_all <- prepare.pheno(pheno_all)
pheno_all <- cbind(pheno_all,
               binary.from.categorical(pheno_all$FCbox,paste0("FCbox",1:4)),
               binary.from.categorical(pheno_all$PPIbox,paste0("PPIbox",1:5)),
               binary.from.categorical(pheno_all$methcage,
                                       paste0("methcage",1:12)),
               binary.from.categorical(pheno_all$round,paste0("SW",1:25)))


# LOAD GENOTYPE DATA
# ------------------
# Load the "mean genotypes"; i.e., the the mean alternative allele
# counts.
cat("Loading genotype data.\n")
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


X_all[1:30,1:10]

map[1:10,]

pheno_all[1:10,1:20]

analysis_selection <- analyses
# min_var = 1
# max_var = 100

# Example for all variants
min_var = 1
max_var = dim(X_all)[2]

# Example to select only one analysis
# analysis_selection <- analyses["BMD"]

# Example to select only one chromosome
# min_var = min(which(map["chr"]==11))
# max_var = max(which(map["chr"]==11))

# Examples from sub-graphs (b) and (c) at https://www.nature.com/articles/ng.3609/figures/3 
# analysis_selection <- analyses["pp12PPIavg"]
# min_var = min(which(map["chr"]==11))
# max_var = max(which(map["chr"]==11))

# analysis_selection <- analyses["pp12PPIavg"]
# min_var = min(which(map["chr"]==7))
# max_var = max(which(map["chr"]==7))

# analysis_selection <- analyses["testis"]
# min_var = min(which(map["chr"]==13))
# max_var = max(which(map["chr"]==13))

# analysis_selection <- analyses["testis"]
# min_var = min(which(map["chr"]==5))
# max_var = max(which(map["chr"]==5))

chromosomes <- unique(map[min_var:max_var,"chr"])

library(qqman)
library(data.table)

for(analysis in analysis_selection) { 

    ##################################
    # Cleanup data
    cat("Loading experiment\n")
    phenotype  <- analysis$pheno
    covariates <- analysis$cov
    outliers   <- analysis$outliers
    
    cat("Removing outliers from experiment\n")
    pheno <- copy(pheno_all)
    if (!is.null(outliers))
      pheno_all <- remove.outliers(pheno,phenotype,covariates,outliers)

    
    # Only analyze samples (i.e. rows of the genotype and phenotype
    # matrices) for which the phenotype and all the covariates are
    # observed.
    cat("Filtering experiment samples with missing data\n")
    pheno <- pheno[which(none.missing.row(pheno[c(phenotype,covariates)])),]

    # Align the phenotypes and genotypes
    cat("Aligning phenotypes and genotypes\n")
    ids   <- intersect(pheno_all$id,rownames(X_all))
    pheno <- pheno_all[match(ids,pheno_all$id),]
    X     <- X_all[match(ids,rownames(X_all)),]

    ###################################
    # Compute using gemma
    # MAP QTLs 
    # Note: Geno is called X, geno and pheno are paired by ID, 
    # pheno has a column called ID, ID in geno is the row number/name.
    X_csv     <- paste(experiment_dir, "/out/X_", analysis$pheno, ".csv",sep="")
    pheno_csv <- paste(experiment_dir, "/out/pheno_", analysis$pheno, ".csv",sep="")
    ids_csv   <- paste(experiment_dir, "/out/ids_", analysis$pheno, ".csv",sep="")
    
    write.csv(X, X_csv, row.names = TRUE)
    write.csv(pheno, pheno_csv, row.names = FALSE)
    write.csv(ids, ids_csv, row.names = FALSE)

    ge_out_dat <- paste(experiment_dir, "/out/ge_", analysis$pheno, "_", min_var, "_", max_var, ".dat", sep="")
    ge_out_csv <- paste(experiment_dir, "/out/ge_", analysis$pheno, "_", min_var, "_", max_var, ".csv", sep="")
    ge_out_png <- paste(experiment_dir, "/out/ge_", analysis$pheno, "_", min_var, "_", max_var, ".png", sep="")
    
    if (!file.exists(ge_out_csv) | !file.exists(ge_out_png)) {
        # Calculate p-values using GEMMA.
        gwscan.gemma <- run.gemma(phenotype,covariates,pheno,X,map,
                                  gemmadir,gemma.exe,chromosomes)

        # Save results to file.
        cat("Saving results to file.\n")
        save(list = c("analysis","gwscan.gemma"),file = ge_out_dat)
        
        named_gws <- gwscan.gemma
        named_gws$snp = rownames(named_gws)
        named_gws$p = 10 ^ (-named_gws$log10p)
        
        png(filename = ge_out_png,
            width = 600, height = 600, units = "px", pointsize = 12,
             bg = "white",  res = NA,
             type = c("cairo", "cairo-png", "Xlib", "quartz"))
        manhattan(named_gws, chr="chr", bp="pos", snp='snp', p="p", annotatePval = 0.1, annotateTop = TRUE )
        dev.off()
        
        write.csv(data.table(named_gws)[order(rank(p)),], ge_out_csv)
        
    }
    
    ###################################
    # Compute using linear model
    lm_out_csv <- paste(experiment_dir, "/out/lm_", analysis$pheno, "_", min_var, "_", max_var, ".csv",sep="")
    lm_out_png <- paste(experiment_dir, "/out/lm_", analysis$pheno, "_", min_var, "_", max_var, ".png",sep="")
    
    if(!file.exists(lm_out_csv) | ! file.exists(lm_out_png)) {
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
        png(filename = lm_out_png,
            width = 600, height = 600, units = "px", pointsize = 12,
             bg = "white",  res = NA,
             type = c("cairo", "cairo-png", "Xlib", "quartz"))
        manhattan(dt[min_var:max_var], chr="chr", bp="pos", snp='snp', p="p", annotatePval = 0.1, annotateTop = TRUE)
        write.csv(dt[order(rank(p)),][1:(max_var-min_var)], lm_out_csv)
        dev.off()
        print(paste("Manhattan plot output to: ", lm_out_png, ", (sorted) values saved in: ", lm_out_csv, sep=""))

        ## Print also on screen
        manhattan(dt[min_var:max_var], chr="chr", bp="pos", snp='snp', p="p", annotatePval = 0.1, annotateTop = TRUE)
    }
}


