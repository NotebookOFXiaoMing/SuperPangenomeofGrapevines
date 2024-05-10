# bwa index ../17.callSVusingUpdateRef/VHP.hap1.fna -p ref.index/vhpHap1

SAMPLES, = glob_wildcards("01.clean.fq/{sample}_1.fq.gz")

print("total sample size: ",len(SAMPLES))

rule all:
    input:
        expand("02.sorted.bam/{sample}.sorted.bam",sample=SAMPLES)

rule bwa:
    input:
        r1 = "01.clean.fq/{sample}_1.fq.gz",
        r2 = "01.clean.fq/{sample}_2.fq.gz",
        ref = "../17.callSVusingUpdateRef/VHP.hap1.fna"
    output:
        "02.sam/{sample}.sam"
    threads:
        16
    resources:
        mem_mb = 96000
    params:
        "ref.index/vhpHap1"
    shell:
        """
        bwa mem -t {threads} {params} {input.r1} {input.r2} -o {output}
        """

rule samtools_sort:
    input:
        rules.bwa.output
    output:
        "02.sorted.bam/{sample}.sorted.bam"
    threads:
        16
    resources:
        mem_mb = 64000
    shell:
        """
        samtools sort -@ {threads} -O bam -o {output} {input}
        """