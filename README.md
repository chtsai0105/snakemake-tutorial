## Install the workflow manage system [**snakemake**](https://snakemake.readthedocs.io/en/stable/index.html)
If you're on UCR hpcc you can simply load snakemake from enviornment module.
```
module load snakemake
```

If you're NOT on UCR hpcc and you don't have snakemake in the enviornment module, please follow the below steps to create a snakemake conda environment.
1. First, create an environment named **snakemake**.

    ```
    conda create -n snakemake
    ```

    After create the environment, activate it by:
    
    ```
    conda activate snakemake
    ```

2. Install the package **mamba**, which is a faster version of **conda**. 

    ```
    conda install -c conda-forge mamba
    ```
    
    After **mamba** being installed, you can later switch from `conda install [package]` to `mamba install [package]` to speed up the package installation.

3. Next, install the package **snakemake** through **mamba**.
    
    ```
    mamba install snakemake
    ```
    
4. Then you can execute the `snakemake --help` to show the snakemake helping page. Snakemake is a Python based language and execution environment for GNU Make-like workflows. The workflow is defined by specifying rules in a `snakefile`. Rules further specify how to create sets of output files from sets of input files as well as the parameters and the command. Snakemake automatically determines the dependencies between the rules by matching file names. Please refer to [snakemake tutorial page](https://snakemake.readthedocs.io/en/stable/tutorial/basics.html) to see how to define the rules.

## Clone the metagenome workflow

Clone the repo to your computer.

Clone by the following command if you're using public key for github connection.

```
git clone --recurse-submodules git@github.com:chtsai0105/metagenome-snakemake.git
```

Or clone by https link.

```
git clone --recurse-submodules https://github.com/chtsai0105/metagenome-snakemake.git
```

Otherwise, clone the submodules as a second step by:
```
cd metagenome-snakemake
git submodule update --init
```

Next, go to the directory by `cd metagenome-snakemake`. It should contains the following files:

File                    |Description
------------------------|---------------------------------
`snakefile`             |Define the rules for the workflow.
`metadata.csv`          |The metadata for samples. Define the names of the samples and the fastq files.
`run_snakemake.bash`    |The bash script for running the workflow.
`data/`                 |The folder for the data and the workflow outputs.
`envs/`                 |The folder that contains the yaml config for conda environments.
`slurm/`                |The folder that contains the slurm profile for stajichlab partition@UCR hpcc.

## Define the samples

You should properly defined your metadata, which is recorded in the `metadata.csv`, before running the workflow. There are 2 columns in this csv table - **sample** and **accession**. The column **sample** defined the sample name. You can change it to names which are more distinguishable instead of accession numbers. The column **accession** defined the fastq file GEO accession number you placed in the folder `data/fastq`. Please make sure they are identical to the fastq files you have otherwise the workflow may have trouble to input the files. Please also confirm that the names in each column are unique.

## Setup the slurm profile (Optional)

If you want to run the workflow on the cluster, you have to setup the profile first. Since the `slurm` profile have been cloned as a submodule, you can simply symlink the **slurm** folder to your snakemake config folder by:

```
ln -s `realpath slurm` $HOME/.config/snakemake
```

Please refer to [snakemake profile for slurm](https://github.com/chtsai0105/snakemake_profile-slurm) for additional info.

## Run the workflow

After compiling the template and setup the paramters, the next step is to run the workflow.

Snakemake provide a dry-run feature to examine the workflow before truly running it. You should always test the workflow beforehand to make sure it execute as expected by the following command:

```
snakemake -np
```

After confirming all the steps. You can run the workflow by executing the script `run_snakemake.bash` or the following command:

```
snakemake -p --profile slurm --use-envmodules --use-conda --jobs 2 --max-threads 8
```

Currently I constrain the cpu usage of all parallel jobs to 8, you may change it as you like.
