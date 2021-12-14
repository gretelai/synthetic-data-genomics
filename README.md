# Synthetic Data Genomics
The code in this repository uses Gretel.ai's synthetic data APIs to create synthetic (artificial) versions of real world mouse genotype and connected phenotype datasets. We then measure the accuracy of our synthetic data by replicating the results of a Genome Wide Association Study (GWAS) on the real world genotypes and phenotypes for 1,220 mice from this paper: https://doi.org/10.1038/ng.3609. 

## Installation

Requirements:
* Conda package manager
* Ubuntu 18.04 recommended
* NVidia T4 or faster GPU
* Gretel.ai API key (https://console.gretel.cloud)


Install the Conda package manager:

```
conda create --name genomics python=3.9
conda activate genomics
conda install jupyter
```

## Recreate the original paper experiments
Follow the steps in `EXPERIMENTS.md` to download the experiment datasets and recreate the results from the paper using real world data.

## Synthesize genome and phenome data, run experiments
Next, create synthetic versions of the mouse phenome and genome datasets from the original experiments.
1. `synthetics/01_create_phenome_training_data.ipynb` creates the genome training set and filter irrelevant fields.
2. `synthetics/02_create_synthetic_mouse_phenomes.ipynb` trains a synthetic model on the mouse phenome set.
3. `synthetics/03_build_genome_training_set.ipynb` creates a genome dataset based on abBMD SNPs 
4. `synthetics/04_create_synthetic_mouse_genomes.ipynb` trains a synthetic model on the mouse genome set, runs GWAS analysis and compares to original results

## Additional resources
* `research_paper_code/notebooks/05_compare_associations.ipynb` compute precision, recall and F1 scores for the final synthetic data
* `synthetics/Optional_tune_synthetic_training_params` optionally use Optuna to tune synthetic training parameters.
* `research_paper_code/notebooks/Manhattan plot.ipynb` compute Manhattan plots for both the original and synthetic genome/phenome gwas p-values



