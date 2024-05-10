SAMPLES, = glob_wildcards("/data/wangj/myan/grape/00.illumina.reads/{sample}_R1.fq.gz")

print("Total sample size: ",len(SAMPLES))


rule all:
    input:
        expand("01.clean.fq/{sample}_1.fq.gz",sample=SAMPLES)

rule fastp:
    input:
        r1 = "/data/wangj/myan/grape/00.illumina.reads/{sample}_R1.fq.gz",
        r2 = "/data/wangj/myan/grape/00.illumina.reads/{sample}_R2.fq.gz"
    output:
        r1 = "01.clean.fq/{sample}_1.fq.gz",
        r2 = "01.clean.fq/{sample}_2.fq.gz",
        html = "01.fastp.report/{sample}.html",
        json = "01.fastp.report/{sample}.json"
    threads:
        16
    resources:
        mem = 24000
    shell:
        """
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -w {threads} -h {output.html} -j {output.json}
        """