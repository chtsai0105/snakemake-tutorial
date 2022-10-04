#!/bin/bash

snakemake -p --use-envmodules --use-conda --profile ./snakemake_profile-slurm --jobs 2
