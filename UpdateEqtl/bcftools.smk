SAMPLES, = glob_wildcards("03.sorted.md.bam/{sample}.sorted.md.bam")

print("Total sample size: ",len(SAMPLES))


rule all:
    input:
        expand("05.vcf/{sample}.vcf",sample=SAMPLES)


rule bcftools_mpileup:
    input:
        bam = "03.sorted.md.bam/{sample}.sorted.md.bam",
        ref = "/data/wangj/myan/grape/17.callSVusingUpdateRef/VHP.hap1.fna"
    output:
        "04.bcf/{sample}.bcf"
    threads:
        26
    resources:
        mem_mb = 96000
    shell:
        """
        bcftools mpileup -q 20 -Q 20 -C 50 --fasta-ref {input.ref} {input.bam} -O u -o {output} --threads {threads}
        """

rule bcftools_call:
    input:
        rules.bcftools_mpileup.output
    output:
        "05.vcf/{sample}.vcf"
    threads:
        26
    resources:
        mem_mb = 64000
    shell:
        """
        bcftools call -c -V indels -v {input} -O v -o {output} --threads {threads}
        """