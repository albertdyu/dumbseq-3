configfile: "config.yaml"

rule all:
	input:
	    "results/NormalizedCounts.csv",
        "results/NormOverOne.csv",
		expand("stringtie/mergedcounts/{sample}.gtf",sample=config["samples"]),
		expand("bws/TSS/{sample}.pos.TSS.bw", sample=config["samples"]),
        expand("bws/TSS/{sample}.neg.TSS.bw", sample=config["samples"]),
        expand("bws/paired/{sample}.pos.paired.bw", sample=config["samples"]),
        expand("bws/paired/{sample}.neg.paired.bw", sample=config["samples"])

rule bam_to_bigwig_dedup:
    input:
        bam="dedup/{sample}.dedup.bam",
        norm_factors="results/norm_factors.txt"
    output:
        pos_bw="bws/paired/{sample}.pos.paired.bw",
        neg_bw="bws/paired/{sample}.neg.paired.bw"
    threads: 8
    shell:
        """
        set -o pipefail
    	matched_line=$(awk -F'\t' -v sample="dedup.{wildcards.sample}.dedup.bam" '($1 == sample) {{print}}' {input.norm_factors})
    	echo "Matched line: $matched_line"
    	scaleFactor=$(awk -F'\t' -v sample="dedup.{wildcards.sample}.dedup.bam" '($1 == sample) {{print $2}}' {input.norm_factors})
    	echo "scaleFactor: $scaleFactor"
        bamCoverage --bam {input.bam} -bs 1 --outFileFormat bigwig --outFileName {output.pos_bw} --filterRNAstrand forward --scaleFactor $scaleFactor -p {threads}
        bamCoverage --bam {input.bam} -bs 1 --outFileFormat bigwig --outFileName {output.neg_bw} --filterRNAstrand reverse --scaleFactor $scaleFactor -p {threads}
        """
		
rule bam_to_bigwig_TSS:
    input:
        bam="TSS/{sample}.TSS_sorted.bam",
        norm_factors="results/norm_factors.txt"
    output:
        pos_bw="bws/TSS/{sample}.pos.TSS.bw",
        neg_bw="bws/TSS/{sample}.neg.TSS.bw"
    shell:
        """
        set -o pipefail
    	scaleFactor=$(awk -F'\t' -v sample="dedup.{wildcards.sample}.dedup.bam" '($1 == sample) {{print $2}}' {input.norm_factors})
        bamCoverage --bam {input.bam} -bs 1 --outFileFormat bigwig --outFileName {output.pos_bw} --filterRNAstrand forward --scaleFactor $scaleFactor -p 8
        bamCoverage --bam {input.bam} -bs 1 --outFileFormat bigwig --outFileName {output.neg_bw} --filterRNAstrand reverse --scaleFactor $scaleFactor -p 8
        """

rule deseq2_normalization:
    input:
        counts="results/counts.featureCounts"
    output:
        "results/norm_factors.txt",
        "results/NormalizedCounts.csv",
        "results/NormOverOne.csv"
    shell:
        "Rscript scripts/deseq2_normalization.R {input.counts} {output}"

rule feature_counts:
    input:
        samples=expand("dedup/{sample}.dedup.bam", sample=config["samples"]), # list of sam or bam files
        annotation=config["annotation"],
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    output:
        multiext(
            "results/{sample}",
            ".featureCounts",
            ".featureCounts.summary",
            ".featureCounts.jcounts",
        ),
    threads: 2
    params:
        strand=0,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        r_path="",  # implicitly sets the --Rpath flag
        extra="-O --fracOverlap 0.2 -J -p",
    log:
        "logs/{sample}.log",
    wrapper:
        "v1.27.0-17-g30c2a724/bio/subread/featurecounts"
        
rule stringtie_merged_count:
	input:
		bam="dedup/{sample}.dedup.bam",
		gtf="stringtie/stringtiemerged.gtf"
	output:
		"stringtie/mergedcounts/{sample}.gtf"
	shell:
		"stringtie -e -B -p 8 -G {input.gtf} -o {output} {input.bam}"
        
rule stringtie_merge:
	input:
		gtf=expand("stringtie/{sample}.gtf",sample=config["samples"]),
		annotation=config["annotation"]
	output:
		"stringtie/stringtiemerged.gtf"
	shell:
		"stringtie --merge -i -p 8 -G {input.annotation} -o {output} {input.gtf}" 

rule stringtie:
	input:
		bam="dedup/{sample}.dedup.bam",
		annotation=config["annotation"]		
	output:
		"stringtie/{sample}.gtf"
	shell:
		"stringtie -p 8 -o {output} -G {input.annotation} --fr -m 100 {input.bam}"

rule getTSS:
	input:
		bam="clipped/{sample}.clipped.bam",
		bai="clipped/{sample}.clipped.bam.bai"
	output:
		"TSS/{sample}.TSS_sorted.bam",
		"TSS/{sample}.TSS_sorted.bam.bai"
	shell:
		"python scripts/get_TSS_bam.py -d ./TSS/ -s {wildcards.sample}.TSS -f {input.bam}"

rule indexclipped:
    input:
        "clipped/{sample}.clipped.bam"
    output:
        temp("clipped/{sample}.clipped.bam.bai")
    log:
        "logs/samtools_index/{sample}.dedup.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        1     # This value - 1 will be sent to -@
    wrapper:
        "0.79.0/bio/samtools/index"

rule removeclipping:
	input:
		bam="dedup/{sample}.dedup.bam",
		bai="dedup/{sample}.dedup.bam.bai"
	output:
		temp("clipped/{sample}.clipped.bam")
	shell:
		"python scripts/removeclipping.py {input.bam} {output}"
		
rule samtools_index_dedup:
    input:
        "dedup/{sample}.dedup.bam"
    output:
        "dedup/{sample}.dedup.bam.bai"
    log:
        "logs/samtools_index/{sample}.dedup.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        1     # This value - 1 will be sent to -@
    wrapper:
        "0.79.0/bio/samtools/index"
        
rule umitools_dedup:
	input:
		bam="sorted/{sample}.sorted.bam",
		bai="sorted/{sample}.sorted.bam.bai"
	output:
		"dedup/{sample}.dedup.bam"
	shell:
		"umi_tools dedup -I {input.bam} --paired -S {output}"

rule samtools_index:
    input:
        "sorted/{sample}.sorted.bam"
    output:
        "sorted/{sample}.sorted.bam.bai"
    log:
        "logs/samtools_index/{sample}.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@
    wrapper:
        "0.79.0/bio/samtools/index"

rule samtools_sort:
    input:
        "aligned/{sample}.bam"
    output:
        "sorted/{sample}.sorted.bam"
    params:
        extra = "-m 4G",
        tmp_dir = "/tmp/"
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@.
    wrapper:
        "0.79.0/bio/samtools/sort"

rule samtools_view:
	input:
		 "aligned/{sample}/Aligned.out.sam"
	output:
		temp("aligned/{sample}.bam")
	threads: 4
	shell:
		"samtools view -@ {threads} -bhS {input} > {output}"
		
rule star_pe_multi:
    input:
        fq1 = "umiextracted_reads/{sample}_R1.UMIextract.fq.gz",
        fq2 = "umiextracted_reads/{sample}_R2.UMIextract.fq.gz",
        idx = config["index"],
    output:
        aln=temp("aligned/{sample}/Aligned.out.sam"),
        log="logs/pe/{sample}/Log.out",
        sj="aligned/{sample}/SJ.out.tab"
    log:
        "logs/star/pe/{sample}.log"
    params:
        # optional parameters
        extra="--limitSjdbInsertNsj 2000000 --outFilterMultimapNmax 1"
    threads: 8
    wrapper:
        "v1.25.0/bio/star/align"
        
rule umiextract:
    input:
        read1="trimmed/{sample}_1.trimmed.fastq.gz",
        read2="trimmed/{sample}_2.trimmed.fastq.gz",
    output:
        out1="umiextracted_reads/{sample}_R1.UMIextract.fq.gz",
        out2="umiextracted_reads/{sample}_R2.UMIextract.fq.gz",
    shell:
        "umi_tools extract -I {input.read1} --extract-method=regex --bc-pattern='^(?P<umi_1>.{{8}})' --read2-in={input.read2} --stdout={output.out1} --read2-out={output.out2}"

rule cutadapt:
    input:
        fq1=lambda wildcards: "{}_R1_001.fastq.gz".format(config["samples"][wildcards.sample], wildcards.sample),
        fq2=lambda wildcards: "{}_R2_001.fastq.gz".format(config["samples"][wildcards.sample], wildcards.sample),
    output:
        trimmed1=temp("trimmed/{sample}_1.trimmed.fastq.gz"),
        trimmed2=temp("trimmed/{sample}_2.trimmed.fastq.gz"),
    threads: 4
    shell:
        "cutadapt -g ATTGCGCAATG --discard-untrimmed --cores={threads} -m 20 -o {output.trimmed1} -p {output.trimmed2} {input.fq1} {input.fq2}"