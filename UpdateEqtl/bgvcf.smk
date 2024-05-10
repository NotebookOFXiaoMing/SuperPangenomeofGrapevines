#SAMPLES, = glob_wildcards("05.vcf/{sample}.vcf")

SAMPLES = []

with open("113samples.id",'r') as fr:
    for line in fr:
        SAMPLES.append(line.strip())

print("Total sample size: ",len(SAMPLES))

rule all:
    input:
        expand("07.changeSampleID.vcf/{sample}.vcf.gz.tbi",sample=SAMPLES),
        "08.merged.vcf/merged.all.vcf.gz.tbi"


rule header:
    input:
        "05.vcf/{sample}.vcf"
    output:
        "06.vcf.header/{sample}.header"
    threads:
        8
    resources:
        mem_mb = 16000
    shell:
        """
        grep "^#" {input} > {output}
        """

rule changeSampleID:
    input:
        rules.header.output
    output:
        "06.vcf.header/{sample}.header.changeSampleID"
    threads:
        4
    resources:
        mem_mb = 8000
    run:
        fr = open(input[0],'r')
        fw = open(output[0],'w')
        for line in fr:
            if line.startswith("#CHROM"):
                fw.write("%s\n"%(line.strip().replace("03.sorted.md.bam/","").replace(".sorted.md.bam","")))
            else:
                fw.write(line)

        fr.close()
        fw.close()

rule body:
    input:
        "05.vcf/{sample}.vcf"
    output:
        "06.vcf.body/{sample}.body"
    threads:
        8
    resources:
        mem_mb = 16000
    shell:
        """
        grep -v "^#" {input} > {output}
        """

rule catHeaderBody:
    input:
        header = rules.changeSampleID.output,
        body = rules.body.output
    output:
        "07.changeSampleID.vcf/{sample}.vcf"
    threads:
        2
    resources:
        mem_mb = 4000
    shell:
        """
        cat {input.header} {input.body} > {output}
        """

rule bgzip:
    input:
        "07.changeSampleID.vcf/{sample}.vcf"
    output:
        "07.changeSampleID.vcf/{sample}.vcf.gz"
    threads:
        8
    resources:
        mem_mb = 16000
    shell:
        """
        bgzip --threads {threads} -c {input} > {output}
        """

rule tabix:
    input:
        rules.bgzip.output
    output:
        "07.changeSampleID.vcf/{sample}.vcf.gz.tbi"
    threads:
        8
    resources:
        mem_mb = 16000
    shell:
        """
        tabix {input}
        """


rule bcftools_merge:
    input:
        vcf = expand(rules.bgzip.output,sample=SAMPLES),
        tbi = expand(rules.tabix.output,sample=SAMPLES)
    output:
        "08.merged.vcf/merged.all.vcf"
    threads:
        96
    resources:
        mem_mb = 230000
    shell:
        """
        /data/wangj/myan/biotools/bcftools-1.17/bcftools merge --missing-to-ref {input.vcf} -o {output} --threads {threads}
        """

rule bgzip_merged:
    input:
        rules.bcftools_merge.output
    output:
        "08.merged.vcf/merged.all.vcf.gz"
    threads:
        8
    resources:
        mem_mb = 16000
    shell:
        """
        bgzip -c {input} > {output}
        """

rule tabix_merged:
    input:
        rules.bgzip_merged.output
    output:
        "08.merged.vcf/merged.all.vcf.gz.tbi"
    threads:
        8
    resources:
        mem_mb = 16000
    shell:
        """
        tabix {input}
        """