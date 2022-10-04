import pandas as pd

## Use the pandas to load the metadata and generate a sample list.
sample_df = pd.read_csv('metadata.csv')
sample_list = sample_df['accession'].tolist()


## By default snakemake executes the first rule in the snakefile. This is a pseudo-rule which is used to define build-targets.
rule all:
    input:
        expand("data/fastqc/pre/{sample}_1_fastqc.zip", sample=sample_list),
        expand("data/fastqc/pre/{sample}_2_fastqc.zip", sample=sample_list),
        expand("data/fastqc/post/{sample}_1_fastqc.zip", sample=sample_list),
        expand("data/fastqc/post/{sample}_2_fastqc.zip", sample=sample_list),
        expand("data/bwa/{sample}.bam", sample=sample_list),
        "data/variance_calling/SARS-CoV2.stat",
        "data/variance_calling/plots"

        ## The above is a more flexible version psuedo rule to replace the hard-coded version below
        # "data/fastqc/pre/SRR17309642_1_fastqc.zip", "data/fastqc/pre/SRR17309642_2_fastqc.zip",
        # "data/fastqc/pre/SRR17309643_1_fastqc.zip", "data/fastqc/pre/SRR17309643_2_fastqc.zip",
        # "data/fastqc/post/SRR17309642_1_fastqc.zip", "data/fastqc/post/SRR17309642_2_fastqc.zip",
        # "data/fastqc/post/SRR17309643_1_fastqc.zip", "data/fastqc/post/SRR17309643_2_fastqc.zip",
        # "data/bwa/SRR17309642.bam", "data/bwa/SRR17309643.bam",
        # "data/variance_calling/SARS-CoV2.stat",
        # "data/variance_calling/plots"


## This is the first "real" rule to run some commands.
## Inside a rule, you should define:
## input: The input files.
## output: The output files. Should be explicitly defined since snakemake will do a check to see whether all the outputs are existed after the execution of a rule.
## params: Optional. You can use for some other parameters. This is especially useful for programs that use output dir rather than explicit output files.
## Threads: How many cpu to use in this rule.
## conda: The yaml files that used to create the conda environment. (activated by flag --use-conda)
## envmodules: Specify the environment module. (activated by flag --use-envmodules)
## shell: The commands or programs to be run. The input, output, params, threads can be specified inside the curly brackets. (e.g. {input})
rule fastqc:
    input:
        "data/fastq/{sample}_1.fastq",      # Make sure to have cammas (,) between each item.
        "data/fastq/{sample}_2.fastq"
    output:
        "data/fastqc/pre/{sample}_1_fastqc.zip",
        "data/fastqc/pre/{sample}_1_fastqc.html",
        "data/fastqc/pre/{sample}_2_fastqc.zip",
        "data/fastqc/pre/{sample}_2_fastqc.html"
    params:
        "data/fastqc/pre"
    threads: 4
    conda:
        "envs/preprocess.yaml"
    envmodules:
        "fastqc/0.11.9"
    shell:
        """
        fastqc -t {threads} -o {params} {input}
        """


rule trimmomatic:
    ## You can use a dict-like objects instead of a list-like objects for input and output rules.
    input:
        R1 = "data/fastq/{sample}_1.fastq",
        R2 = "data/fastq/{sample}_2.fastq"
    output:
        R1_paired = "data/trimmed_fastq/{sample}_1.fastq",
        R1_unpaired = "data/trimmed_fastq/{sample}_1_unpaired.fastq",
        R2_paired = "data/trimmed_fastq/{sample}_2.fastq",
        R2_unpaired = "data/trimmed_fastq/{sample}_2_unpaired.fastq"
    threads: 4
    conda:
        "envs/preprocess.yaml"
    envmodules:
        "trimmomatic/0.39"
    shell:
        ## You can specify the the output by its keys
        """
        trimmomatic PE -threads {threads} {input.R1} {input.R2} {output.R1_paired} {output.R1_unpaired} {output.R2_paired} {output.R2_unpaired} \
        ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15
        """

## Inherit from previously defined rules
## We want to run fastqc again against the trimmed fastq. We can simply reuse the already-defined rule "fastqc"
## and overwrite all the path-dependent parameters. (which is input, output and params)
## More details in https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#rule-inheritance
use rule fastqc as fastqc_post with:
    input:
        "data/trimmed_fastq/{sample}_1.fastq",
        "data/trimmed_fastq/{sample}_2.fastq"
    output:
        "data/fastqc/post/{sample}_1_fastqc.zip",
        "data/fastqc/post/{sample}_1_fastqc.html",
        "data/fastqc/post/{sample}_2_fastqc.zip",
        "data/fastqc/post/{sample}_2_fastqc.html"
    params:
        "data/fastqc/post"


## Use the expand function to make the code more simple.
## More details in https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#the-expand-function
rule bwa_index:
    input:
        "data/ref/GCF_009858895.2_ASM985889v3_genomic.fna.gz"
    output:
        # "data/ref/bwa_idx.0123",
        # "data/ref/bwa_idx.amb",
        # "data/ref/bwa_idx.ann",
        # "data/ref/bwa_idx.bwt.2bit.64",
        # "data/ref/bwa_idx.pac"
        expand("data/ref/bwa_idx.{ext}", ext=["0123", "amb", "ann", "bwt.2bit.64", "pac"])        ## This expression is equivalent to the above 5 lines
    params:
        "data/ref/bwa_idx"
    conda:
        "envs/alignment.yaml"
    envmodules:
        "bwa-mem2/2.2.1"
    shell:
        """
        bwa-mem2 index -p {params} {input}
        """


rule bwa:
    input:
        # expand("data/ref/bwa_idx.{ext}", ext=["0123", "amb", "ann", "bwt.2bit.64", "pac"])
        # "data/trimmed_fastq/{sample}_1.fastq",
        # "data/trimmed_fastq/{sample}_2.fastq"
        idx = rules.bwa_index.output,                   ## Using rule dependencies expression to replace the above 3 lines
        R1 = rules.trimmomatic.output.R1_paired,        ## More details in https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#rule-dependencies
        R2 = rules.trimmomatic.output.R2_paired
    output:
        "data/bwa/{sample}.sam"
    params:
        "data/ref/bwa_idx"
    threads: 4
    conda:
        "envs/alignment.yaml"
    envmodules:
        "bwa-mem2/2.2.1"
    shell:
        """
        bwa-mem2 mem -t {threads} {params} {input.R1} {input.R2} > {output}
        """


rule samtools_fixmate:
    input:
        rules.bwa.output
    output:
        fixmate = temp("data/bwa/{sample}_fixmate.bam"),    ## Define a temporary files, which will be removed after all rules that use it as an input are completed.
        bam = "data/bwa/{sample}.bam"
    threads: 4
    conda:
        "envs/alignment.yaml"
    envmodules:
        "samtools/1.16"
    shell:
        """
        samtools fixmate -O bam {input} {output.fixmate}
        samtools sort --threads {threads} -O BAM -o {output.bam} {output.fixmate}
        samtools index {output.bam}
        """


rule unzip_genome:
    input:
        "data/ref/GCF_009858895.2_ASM985889v3_genomic.fna.gz"
    output:
        temp("data/ref/genome.fna")
    shell:
        """
        gunzip -c {input} > {output}
        """


rule bcftools_pileup:
    input:
        genome = rules.unzip_genome.output,
        bams = expand("data/bwa/{sample}.bam", sample=sample_list)
    output:
        genome_idx = temp("data/ref/genome.fna.fai"),
        vcf = "data/variance_calling/SARS-CoV2.vcf.gz",
        tbi = "data/variance_calling/SARS-CoV2.vcf.gz.tbi"
    conda:
        "envs/variance_calling.yaml"
    envmodules:
        "bcftools/1.16"
    shell:
        """
        bcftools mpileup -Ou -f {input.genome} {input.bams} | bcftools call -vmO z -o {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule bcftools_stats:
    input:
        vcf = rules.bcftools_pileup.output.vcf,
        genome = rules.unzip_genome.output,
        genome_idx = rules.bcftools_pileup.output
    output:
        stats = "data/variance_calling/SARS-CoV2.stat"
    threads: 8
    conda:
        "envs/variance_calling.yaml"
    envmodules:
        "bcftools/1.16"
    shell:
        """
        bcftools stats --threads {threads} -F {input.genome} -s - {input.vcf} > {output.stats}
        """


rule plot:
    input:
        rules.bcftools_stats.output.stats
    output:
        directory("data/variance_calling/plots")
    conda:
        "envs/variance_calling.yaml"
    shell:
        """
        plot-vcfstats -p {output} {input}
        """
