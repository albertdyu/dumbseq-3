# dumbseq-3
A snakemake pipeline for Smart-seq3 analysis

This pipeline does UMI deduplication, mapping, and generates normalized BW files with one command. It will also output some stuff for circadian analysis and does transcript assembly with Stringtie. The latter I haven't messed around with much, but maybe you'll find it useful.

Run the command to create the environment:

conda create --name ss3 -c bioconda cutadapt umi_tools samtools stringtie deeptools subread python r-base snakemake

Activate the environment: 

conda activate ss3

Install python dependencies:

pip install pysam

Install R dependencies:

R
BiocManager::install("DESeq2")
Q()

Edit the config file

![CleanShot 2024-02-15 at 12 16 19@2x](https://github.com/albertdyu/dumbseq-3/assets/13254772/3474cd88-4067-423a-8e78-3a6689078e97)
Run snakemake!
nohup snakemake --cores 8 --configfile config.yaml &
![image](https://github.com/albertdyu/dumbseq-3/assets/13254772/08f55860-5b7b-45cf-b29d-5bc686b05bb9)
