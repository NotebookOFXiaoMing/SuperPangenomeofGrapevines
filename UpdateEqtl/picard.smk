SAMPLES, = glob_wildcards("02.sorted.bam/{sample}.sorted.bam")

print("Total sample size: ",len(SAMPLES))

rule all:
    input:
        expand("03.sorted.md.bam/{sample}.sorted.md.bam",sample=SAMPLES)


rule picard:
    input:
        "02.sorted.bam/{sample}.sorted.bam"
    output:
        bam = "03.sorted.md.bam/{sample}.sorted.md.bam",
        txt = "03.sorted.md.bam/{sample}.sorted_dup_metrics.txt"
    threads:
        32
    resources:
        mem_mb = 128000
    shell:
        """
        picard MarkDuplicates I={input} O={output.bam} Remove_Duplicates=true ASSUME_SORTED=true M={output.txt}
        """