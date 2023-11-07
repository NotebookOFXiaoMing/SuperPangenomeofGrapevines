SAMPLES, = glob_wildcards(config['r1'])

print(SAMPLES)
print("Total sample size: ",len(SAMPLES))

rule all:
    input:
        expand("09.vcf/{sample}.vcf",sample=SAMPLES),
        "merged.vcf.gz"


rule vgGiraffeMap:
    input:
        r1 = config['r1'],
        r2 = config['r2'],
        gbz = config['gbz'],
        min = config['min'],
        dist = config['dist']
    output:
        gam = "07.gam/{sample}.gam"
    threads:
        24
    resources:
        mem_mb = 96000
    benchmark:
        "vg.Map.benchmarks/{sample}.txt"
    shell:
        """
        vg giraffe -Z {input.gbz} -m {input.min} -d {input.dist} \
        -f {input.r1} -f {input.r2} --threads {threads} > {output.gam}
        """

rule vgPack:
    input:
        gam = rules.vgGiraffeMap.output.gam,
        gbz = config['gbz']
    output:
        pack = "08.pack/{sample}.pack"
    threads:
        24
    resources:
        mem_mb = 96000
    benchmark:
        "vg.Pack.benchmarks/{sample}.txt"
    shell:
        """
        vg pack -x {input.gbz} -g {input.gam} -Q 5 -s 5 \
        -o {output.pack} -t {threads}
        """

rule vgCall:
    input:
        gbz = config['gbz'],
        pack = rules.vgPack.output.pack
    output:
        vcf = "09.vcf/{sample}.vcf"
    resources:
        mem_mb = 96000
    benchmark:
        "vgPack.benchmarks/{sample}.txt"
    shell:
        """
        vg call {input.gbz} -k {input.pack} -a -t {threads} > {output}
        """

rule rename_vcf:
    input:
        rules.vgCall.output.vcf
    output:
        "10.rename.vcf/{sample}.vcf"
    threads:
        4
    resources:
        mem_mb = 4000
    params:
        "{sample}"
    run:
        fr = open(input[0],'r')
        fw = open(output[0],'w')

        for line in fr:
            if line.startswith("##"):
                fw.write(line)
            elif line.startswith("#CHROM"):
                fw.write("%s\n"%line.strip().replace("SAMPLE",params[0]))
            else:
                fw.write(line)

        fr.close()
        fw.close()


rule bgzip:
    input:
        rules.rename_vcf.output
    output:
        "10.rename.vcf/{sample}.vcf.gz"
    threads:
        8
    resources:
        mem_mb = 16000
    shell:
        """
        bgzip -c {input} > {output} -@ {threads}
        """

rule tabix:
    input:
        rules.bgzip.output
    output:
        "10.rename.vcf/{sample}.vcf.gz.tbi"
    threads:
        8
    resources:
        mem_mb = 16000
    shell:
        """
        tabix {input}
        """

rule bcftools_merged:
    input:
        vcf = expand("10.rename.vcf/{sample}.vcf.gz",sample=SAMPLES),
        tbi = expand("10.rename.vcf/{sample}.vcf.gz.tbi",sample=SAMPLES)
    output:
        "merged.vcf",
        "merged.vcf.gz",
        "merged.vcf.gz.tbi"
    threads:
        24
    resources:
        mem_mb = 64000
    shell:
        """
        bcftools merge {input.vcf} -o {output[0]} --threads {threads}
        bgzip -c {output[0]} > {output[1]} -@ {threads}
        tabix {output[1]}
        """