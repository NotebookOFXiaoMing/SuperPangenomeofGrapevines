SAMPLES, = glob_wildcards("/data/wangj/myan/grape/00.20.samples.hifi.reads/{sample}.hifi.clean.fq.gz")

SAMPLES.remove("V017")

print("Total sample size: ",len(SAMPLES))

rule all:
    input:
        expand("01.hifi.sam/{sample}.sam",sample=SAMPLES),
        expand("02.hifi.bam/{sample}.sorted.bam",sample=SAMPLES),
        expand("04.hifi.filter.vcf/{sample}.filter.vcf",sample=SAMPLES)

rule minimap2:
    input:
        ref = "/data/wangj/myan/grape/17.callSVusingUpdateRef/VHP.hap1.fna",
        hifi = "/data/wangj/myan/grape/00.20.samples.hifi.reads/{sample}.hifi.clean.fq.gz"
    output:
        "01.hifi.sam/{sample}.sam"
    threads:
        32
    resources:
        mem_mb = 64000
    shell:
        """
        minimap2 -t {threads} -ax map-hifi {input.ref} {input.hifi} > {output}
        """

rule samtools:
    input:
        rules.minimap2.output
    output:
        "02.hifi.bam/{sample}.sorted.bam"
    threads:
        32
    resources:
        mem_mb = 64000
    shell:
        """
        samtools view -@ 8 -bh {input} | samtools sort -@ 8 -O BAM -o {output}
        """

rule samtools_index:
    input:
        rules.samtools.output
    output:
        "02.hifi.bam/{sample}.sorted.bam.bai"
    threads:
        32
    resources:
        mem_mb = 64000
    shell:
        """
        samtools index -@ 8 {input}
        """

rule sniffles:
    input:
        bam = rules.samtools.output,
        bai = rules.samtools_index.output
    output:
        "03.hifi.vcf/{sample}.vcf"
    threads:
        16
    resources:
        mem_mb = 64000
    shell:
        """
        sniffles --threads {threads} --minsvlen 50 --input {input.bam} --vcf {output}
        """

rule filterVCF:
    input:
        vcf = rules.sniffles.output,
        ref = "/data/wangj/myan/grape/17.callSVusingUpdateRef/VHP.hap1.fna"
    output:
        "04.hifi.filter.vcf/{sample}.filter.vcf"
    threads:
        16
    resources:
        mem_mb = 32000
    params:
        "{sample}"
    shell:
        """
        python editSnifflesVCF.py {input.ref} {input.vcf} {output} {params}
        """