import glob

SAMPLES = set([sample.split("/")[-1].replace(".txt","") for sample in glob.glob("00.NLRannotator.output/*/*.txt")])

print(SAMPLES)

print("Total sample size: ",len(SAMPLES))


rule all:
    input:
        expand("01.NLRannoResultsBed/{sample}.hap2/{sample}.bed",sample=SAMPLES),
        expand("01.NLRannoResultsBed/{sample}.hap1/{sample}.bed",sample=SAMPLES),
        expand("01.NLRannoResultsBed/{sample}.hap2/{sample}.sorted.bed",sample=SAMPLES),
        expand("01.NLRannoResultsBed/{sample}.hap1/{sample}.sorted.bed",sample=SAMPLES),
        expand("02.candidate.NLR/{sample}.hap1/{sample}.candidate.NLR.id",sample=SAMPLES),
        expand("02.candidate.NLR/{sample}.hap2/{sample}.candidate.NLR.id",sample=SAMPLES)

rule hap1:
    input:
        "00.NLRannotator.output/{sample}.hap1/{sample}.txt"
    output:
        "01.NLRannoResultsBed/{sample}.hap1/{sample}.bed"
    threads:
        2
    resources:
        mem_mb = 4000
    shell:
        """
        cat {input} | awk '{{print $1"\\t"$4"\\t"$5"\\t"$3}}' > {output}
        """


rule hap2:
    input:
        "00.NLRannotator.output/{sample}.hap2/{sample}.txt"
    output:
        "01.NLRannoResultsBed/{sample}.hap2/{sample}.bed"
    threads:
        2
    resources:
        mem_mb = 4000
    shell:
        """
        cat {input} | awk '{{print $1"\\t"$4"\\t"$5"\\t"$3}}' > {output}
        """


rule bedtools_sort_hap1:
    input:
        fna_bed = "01.gff2bed/{sample}.hap1/{sample}.fna.bed",
        gene_bed = rules.hap1.output
    output:
        "01.NLRannoResultsBed/{sample}.hap1/{sample}.sorted.bed"
    threads:
        4
    resources:
        mem_mb = 4000
    shell:
        """
        bedtools sort -g {input.fna_bed} -i {input.gene_bed} > {output}
        """

rule bedtools_sort_hap2:
    input:
        fna_bed = "01.gff2bed/{sample}.hap2/{sample}.fna.bed",
        gene_bed = rules.hap2.output
    output:
        "01.NLRannoResultsBed/{sample}.hap2/{sample}.sorted.bed"
    threads:
        4
    resources:
        mem_mb = 4000
    shell:
        """
        bedtools sort -g {input.fna_bed} -i {input.gene_bed} > {output}
        """

rule bedtools_closest_hap1:
    input:
        a = rules.bedtools_sort_hap1.output,
        b = "01.gff2bed/{sample}.hap1/{sample}.sorted.bed",
        g = "01.gff2bed/{sample}.hap1/{sample}.fna.bed"
    output:
        id = "02.candidate.NLR/{sample}.hap1/{sample}.candidate.NLR.id",
        group = "02.candidate.NLR/{sample}.hap1/{sample}.candidate.NLR.group"
    threads:
        2
    resources:
        mem_mb = 4000
    shell:
        """
        bedtools closest -a {input.a} -b {input.b} -g {input.g} -D a | awk \
        '$9==0 {{print $0}}' | awk '{{print $8}}' | sort -u > {output.id}

        bedtools closest -a {input.a} -b {input.b} -g {input.g} -D a | awk \
        '$9==0 {{print $0}}' | awk '{{print $8"\\t"$4}}' > {output.group}
        """


rule bedtools_closest_hap2:
    input:
        a = rules.bedtools_sort_hap2.output,
        b = "01.gff2bed/{sample}.hap2/{sample}.sorted.bed",
        g = "01.gff2bed/{sample}.hap2/{sample}.fna.bed"
    output:
        id = "02.candidate.NLR/{sample}.hap2/{sample}.candidate.NLR.id",
        group = "02.candidate.NLR/{sample}.hap2/{sample}.candidate.NLR.group"
    threads:
        2
    resources:
        mem_mb = 4000
    shell:
        """
        bedtools closest -a {input.a} -b {input.b} -g {input.g} -D a | awk \
        '$9==0 {{print $0}}' | awk '{{print $8}}' | sort -u > {output.id}

        bedtools closest -a {input.a} -b {input.b} -g {input.g} -D a | awk \
        '$9==0 {{print $0}}' | awk '{{print $8"\\t"$4}}' > {output.group}
        """