# Recreate the initial experiments
We will recreate the experiments from "Genome-wide association study of behavioral, physiological and gene expression traits in outbred CFW mice". The command below will download the datasets used in the mouse genome paper into a local directory called `mice_data_set` for processing.

## Download the original experiment datasets
The shell script below downloads the original experiment datasets and results for the `abBMD` analysis.

```
sh download.sh
```

## Set up a Jupyter notebook with the R Kernel
Next, use the Conda package manager to set up a virtual environment to run the Jupyter notebooks that recreate the original experiments on the datasets.

```
conda create -n r-kernel
conda install r-recommended r-irkernel
conda install jupyter
```

Add the R-kernel spec to Jupyter and install required packages.
```
R -e 'IRkernel::installspec()'
R -e 'install.packages("qtl", repos = "http://cran.us.r-project.org")'
R -e 'install.packages("qqman", repos = "http://cran.us.r-project.org")'
R -e 'install.packages("data.table", repos = "http://cran.us.r-project.org")'
R -e 'install.packages("stringr", repos = "http://cran.us.r-project.org")'
R -e 'install.packages("qqman", repos = "http://cran.us.r-project.org")'
R -e 'install.packages("devtools", repos = "http://cran.us.r-project.org")'
```

Run Jupyter notebook

```
jupyter notebook
```

## Run Map.ipynb
Next, open `./research_paper_code/notebooks/map.ipynb` in Jupyter notebook, and choose Kernel->Run All. This will run through the R-studio code in this repository that recreates the results from the original paper. As data is generated, you will see plots and data files generated in the following formats:

```
(base) redlined@redlined-980:~/GitHub/synthetic-data-genomics/mice_data_set/out$ head lm_plantaris_1_79646.csv
"","snp","chr","pos","p"
"1","rs46110548",13,8871303,9.60141540760935e-06
"2","rs216161522",13,9038407,1.04292648303643e-05
"3","rs48201941",13,5974272,1.06707961350016e-05
"4","rs231489766",13,8988261,1.131883198666e-05
"5","rs259152022",13,9038004,1.19756669281965e-05
"6","rs48351148",13,8940993,1.26057184411879e-05
"7","rs250637251",15,101019128,1.35049345326697e-05
"8","rs222853368",13,6648524,1.86357665221831e-05
"9","rs243591869",13,7296440,1.89613649778877e-05
```

![image](https://user-images.githubusercontent.com/6510818/136842534-32ed43ac-e80f-47b2-9788-2f5f7149d257.png)




