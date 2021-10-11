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
Running this notebook before synthesizing data is optional as the original `abBMD` analysis was downloaded in the steps above. Map.ipynb can be run to optinally recreate the original experiment results. To run the notebook, open `./research_paper_code/notebooks/map.ipynb` in Jupyter notebook, and choose Kernel->Run All. This will run through the R-studio code in this repository that recreates the results from the original paper. As data is generated, you will see plots and data files generated in the following formats:

```
(base) redlined@redlined-980:~/GitHub/synthetic-data-genomics/mice_data_set/out$ head lm_abBMD_1_79646.csv
"","snp","chr","pos","p"
"1","rs29477109",11,95292217,5.05231663641996e-14
"2","rs27071351",11,96114911,7.07418067212828e-14
"3","rs27024162",11,96918116,7.17058199722633e-14
"4","rs49423067",11,96918212,7.19866140655625e-14
"5","rs29470802",11,95263588,8.04984862217419e-14
"6","rs29459746",11,95987376,1.03725122425739e-13
"7","rs50417410",11,97011284,1.04333530152468e-13
"8","rs29473466",11,96920033,1.33866242959213e-13
"9","rs221074340",11,96018255,1.35574083178291e-13
```

![image](https://user-images.githubusercontent.com/6510818/136842534-32ed43ac-e80f-47b2-9788-2f5f7149d257.png)




