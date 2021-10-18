# Synthetic Data Genomics
We will create synthetic versions of both the mouse genotype and connected phenotype datasets and prove their utility by replicating the results from the paper: Genome-wide association study of behavioral, physiological and gene expression traits in outbred CFW mice by Parker, et al.


## Installation

Requirements:
* Conda package manager
* Ubuntu 18.04 recommended


Install the Conda package manager:

```
conda create --name genomics python=3.9
conda activate genomics
conda install jupyter
```

## Download the original experiment datasets
The shell script below downloads the original experiment datasets and results for the `abBMD` analysis.

```
sh download.sh
```

## Recreate the original paper experiments
(Optional) Follow the steps in `EXPERIMENTS.md` to recreate the results from the paper.

## Synthesize genome and phenome data
Next, create synthetic versions of the mouse phenome and genome datasets from the original experiments.
1. `synthetics/01_create_phenome_training_data.ipynb` creates the genome training set and filter irrelevant fields.
2. `synthetics/02_create_synthetic_mouse_phenomes.ipynb` trains a synthetic model on the mouse phenome set.
4. `synthetics/03_build_genome_training_set.ipynb` creates a genome dataset based on abBMD SNPs 
5. `synthetics/04_create_synthetic_mouse_genomes.ipynb` trains a synthetic model on the mouse genome set.

## Run experiments and compare results
To be added.

## Citations
Parker, C., Gopalakrishnan, S., Carbonetto, P. et al. Genome-wide association study of behavioral, physiological and gene expression traits in outbred CFW mice. Nat Genet 48, 919–926 (2016). https://doi.org/10.1038/ng.3609

