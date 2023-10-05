SAMPLES, = glob_wildcards("../09.panGeneFamilies.144.hap/00.all.peps/{sample}_hap1.fa")

print("Total sample size: ",len(SAMPLES))


rule all:
    input:
        expand("02.candidate.NLR.peps/{sample}.hap1/{sample}_hap1_nlr.fa",sample=SAMPLES),
        expand("02.candidate.NLR.peps/{sample}.hap2/{sample}_hap2_nlr.fa",sample=SAMPLES)


rule seqkit_grep_hap1:
    input:
        peps = "00.all.peps/{sample}.hap1.protein.fasta",
        ids = "02.candidate.NLR/{sample}.hap1/{sample}.candidate.NLR.id"
    output:
        "02.candidate.NLR.peps/{sample}.hap1/{sample}_hap1_nlr.fa"
    threads:
        4
    resources:
        mem_mb = 4000
    shell:
        """
        seqkit grep -r -f {input.ids} {input.peps} -o {output}
        """

rule seqkit_grep_hap2:
    input:
        peps = "00.all.peps/{sample}.hap2.protein.fasta",
        ids = "02.candidate.NLR/{sample}.hap2/{sample}.candidate.NLR.id"
    output:
        "02.candidate.NLR.peps/{sample}.hap2/{sample}_hap2_nlr.fa"
    threads:
        4
    resources:
        mem_mb = 4000
    shell:
        """
        seqkit grep -r -f {input.ids} {input.peps} -o {output}
        """