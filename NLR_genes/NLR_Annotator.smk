SAMPLES, = glob_wildcards("/data/wangj/annotation/finalgff_72/{sample}.hap1.gff")

print("Total sample size: ",len(SAMPLES))

rule all:
    input:
        expand("00.NLRannotator.output/{sample}.hap1/{sample}.txt",sample=SAMPLES),
        expand("00.NLRannotator.output/{sample}.hap2/{sample}.txt",sample=SAMPLES)


rule nlr_hap1:
    input:
        fa = "/data/wangj/grape_assemblies/finalfna/{sample}.hap1.fna",
        x = "/data/wangj/myan/biotools/NLR-Annotator-NLR-Annotator-2/mot.txt",
        y = "/data/wangj/myan/biotools/NLR-Annotator-NLR-Annotator-2/store.txt"
    output:
        txt = "00.NLRannotator.output/{sample}.hap1/{sample}.txt",
        a = "00.NLRannotator.output/{sample}.hap1/{sample}.motif.fa",
        f = "00.NLRannotator.output/{sample}.hap1/{sample}.flank2k.fa"
    threads:
        48
    resources:
        mem_mb = 96
    shell:
        """
        time java -jar ~/my_data/myan/biotools/NLR-Annotator-NLR-Annotator-2/NLR-Annotator-v2.1b.jar -i \
        {input.fa} -x {input.x} -y {input.y} -o {output.txt} \
        -a {output.a} -f {input.fa} {output.f} 2000 -t {threads}
        """


rule nlr_hap2:
    input:
        fa = "/data/wangj/grape_assemblies/finalfna/{sample}.hap2.fna",
        x = "/data/wangj/myan/biotools/NLR-Annotator-NLR-Annotator-2/mot.txt",
        y = "/data/wangj/myan/biotools/NLR-Annotator-NLR-Annotator-2/store.txt"
    output:
        txt = "00.NLRannotator.output/{sample}.hap2/{sample}.txt",
        a = "00.NLRannotator.output/{sample}.hap2/{sample}.motif.fa",
        f = "00.NLRannotator.output/{sample}.hap2/{sample}.flank2k.fa"
    threads:
        48
    resources:
        mem_mb = 96
    shell:
        """
        time java -jar ~/my_data/myan/biotools/NLR-Annotator-NLR-Annotator-2/NLR-Annotator-v2.1b.jar -i \
        {input.fa} -x {input.x} -y {input.y} -o {output.txt} \
        -a {output.a} -f {input.fa} {output.f} 2000 -t {threads}
        """