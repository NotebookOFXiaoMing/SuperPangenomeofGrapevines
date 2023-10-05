#print(config["raw_fq"]["r1"])

SAMPLES, = glob_wildcards(config["raw_fq"]["r1"])

#print(SAMPLES)

#print("Total sample size: ",len(SAMPLES))


rule all:
    input:
        #expand("01.clean.reads/{sample}_clean_R2.fq",sample=SAMPLES),
        "reference/index/ref.1.ht2",
        expand("03.expression/{sample}/gene_abund.tsv",sample=SAMPLES),
        "04.counts/trans_counts.csv",
        "04.counts/gene_counts.csv"


rule fastp:
    input:
        r1 = config["raw_fq"]["r1"],
        r2 = config["raw_fq"]["r2"]
    output:
        r1 = "01.clean.reads/{sample}_clean_R1.fq",
        r2 = "01.clean.reads/{sample}_clean_R2.fq",
        json = "01.fastp.report/{sample}.json",
        html = "01.fastp.report/{sample}.html"
    threads:
        8
    resources:
        mem_mb = 24000
    params:
        "-f 10 -t 10 -F 10 -T 10"
    benchmark:
        "00.benchmarks/{sample}.fastp.benchmark.txt"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
        -o {output.r1} -O {output.r2} \
        -j {output.json} -h {output.html} \
        -w {threads} {params}
        """


rule hisat2_build:
    input:
        config["reference"]["fa"]
    output:
        "reference/index/ref.1.ht2"
    threads:
        8
    resources:
        mem_mb = 24000
    params:
        "reference/index/ref"
    shell:
        """
        hisat2-build {input} {params}
        """

rule hisat2_align:
    input:
        r1 = rules.fastp.output.r1,
        r2 = rules.fastp.output.r2,
        index = rules.hisat2_build.output
    output:
        "02.sam/{sample}.sam"
    threads:
        8
    resources:
        mem_mb = 48000
    params:
        "reference/index/ref"
    shell:
        """
        hisat2 -p {threads} --dta -x {params} \
        -1 {input.r1} -2 {input.r2} \
        -S {output}
        """

rule samtools_sort:
    input:
        rules.hisat2_align.output
    output:
        "02.bam/{sample}.sorted.bam"
    threads:
        8
    resources:
        mem_mb = 48000
    shell:
        """
        samtools sort -@ {threads} -o {output} {input}
        """

rule samtools_index:
    input:
        rules.samtools_sort.output
    output:
        "02.bam/{sample}.sorted.bam.bai"
    threads:
        8
    resources:
        mem_mb = 48000
    shell:
        """
        samtools index -@ {threads} {input}
        """
    
rule stringtie:
    input:
        gtf = config["reference"]['gtf'],
        bam = rules.samtools_sort.output,
        bai = rules.samtools_index.output
    output:
        trans_gtf = "03.expression/{sample}/transcripts.gtf",
        gene_abund = "03.expression/{sample}/gene_abund.tsv",
        t_data_ctab = "03.expression/{sample}/t_data.ctab"
    threads:
        8
    resources:
        mem_mb = 48000
    shell:
        """
        stringtie -p {threads} -G {input.gtf} -e -B -o {output.trans_gtf} -A {output.gene_abund} {input.bam}
        """

rule tximport:
    input:
        t_data_ctab = expand(rules.stringtie.output.t_data_ctab,sample=SAMPLES)
    output:
        gene_counts = "04.counts/gene_counts.csv",
        trans_counts = "04.counts/trans_counts.csv"
    threads:
        2
    resources:
        mem_mb = 8000
    script:
        "scripts/r/run_tximport.R"
