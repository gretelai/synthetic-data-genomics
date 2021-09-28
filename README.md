# Synthetic Data Genomics
We will create synthetic versions of both the mouse genotype and connected phenotype datasets and prove their utility by replicating the results from this paper: https://doi.org/10.1038/ng.3609. 

## Installation

Requirements:
* Conda package manager
* Ubuntu 18.04 recommended


Install the Conda package manager:

```
conda create --name genomics python=3.9
conda activate genomics
conda install r-core r-recommended r-irkernel jupyter
```

## Recreate the original paper experiments
See `EXPERIMENTS.md`.

## Synthesize genome and phenome data
1. Run `synthetics/01_create_phenome_training_data.ipynb` to create genome training set and filter irrelevant fields.

## Run experiments and compare results
To be added.

