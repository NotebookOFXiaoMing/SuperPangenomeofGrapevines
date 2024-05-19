SAMPLES, = glob_wildcards("/data/wangj/myan/grape/00.20.samples.hifi.reads/{sample}.hifi.clean.fq.gz")

SAMPLES.remove("V017")

print("Total sample size: ",len(SAMPLES))

rule all:
    input:
        expand("05.vg.autoindex/{sample}/{sample}.giraffe.gbz",sample=SAMPLES),
        expand("08.vcf/{sample}.vcf",sample=SAMPLES)

rule bgzip:
    input:
        "/data/wangj/myan/grape/23.SV.genotyping.benchmarks/04.hifi.filter.vcf/{sample}.filter.vcf"
    output:
        "/data/wangj/myan/grape/23.SV.genotyping.benchmarks/04.hifi.filter.vcf/{sample}.filter.vcf.gz"
    threads:
        8
    resources:
        mem_mb = 24000
    shell:
        """
        bgzip {input} -c > {output}
        tabix {output}
        """

rule vgautoindex:
    input:
        ref = "/data/wangj/myan/grape/17.callSVusingUpdateRef/VHP.hap1.fna",
        vcf = rules.bgzip.output
    output:
        gbz = "05.vg.autoindex/{sample}/{sample}.giraffe.gbz",
        min = "05.vg.autoindex/{sample}/{sample}.min",
        dist = "05.vg.autoindex/{sample}/{sample}.dist"
    threads:
        32
    resources:
        mem_mb = 64000
    params:
        folder = "05.vg.autoindex/{sample}",
        prefix = "{sample}"
    shell:
        """
        cd {params.folder}
        if [ -d "tmpdir" ]; then
            /data/wangj/myan/biotools/vg autoindex --workflow giraffe \
            -r {input.ref} \
            -v {input.vcf} -p {params.prefix} -t {threads} -T tmpdir
        else
            mkdir tmpdir
            /data/wangj/myan/biotools/vg autoindex --workflow giraffe \
            -r {input.ref} \
            -v {input.vcf} -p {params.prefix} -t {threads} -T tmpdir
        fi
        """

rule vggiraffe:
    input:
        ref = "/data/wangj/myan/grape/17.callSVusingUpdateRef/VHP.hap1.fna",
        r1 = "/data/wangj/myan/grape/00.illumina.reads/{sample}_R1.fq.gz",
        r2 = "/data/wangj/myan/grape/00.illumina.reads/{sample}_R2.fq.gz",
        gbz = rules.vgautoindex.output.gbz,
        min = rules.vgautoindex.output.min,
        dist = rules.vgautoindex.output.dist
    output:
        "06.gam/{sample}.gam"
    threads:
        32
    resources:
        mem_mb = 64000
    params:
        "{sample}"
    shell:
        """
        /data/wangj/myan/biotools/vg giraffe \
        -Z {input.gbz} \
        -m {input.min} \
        -d {input.dist} \
        -f {input.r1} 
        -f {input.r2} \
        --threads 24 \
        -N {params} > {output}
        """

rule vgpack:
    input:
        ref = "/data/wangj/myan/grape/17.callSVusingUpdateRef/VHP.hap1.fna",
        r1 = "/data/wangj/myan/grape/00.illumina.reads/{sample}_R1.fq.gz",
        r2 = "/data/wangj/myan/grape/00.illumina.reads/{sample}_R2.fq.gz",
        gam = rules.vggiraffe.output,
        gbz = rules.vgautoindex.output.gbz,
        min = rules.vgautoindex.output.min,
        dist = rules.vgautoindex.output.dist
    output:
        "07.pcak/{sample}.pack"
    threads:
        32
    resources:
        mem_mb = 64000
    params:
        "{sample}"
    shell:
        """
        /data/wangj/myan/biotools/vg pack \
        -x {input.gbz} \
        -g {input.gam} \
        -Q 5 -s 5 \
        -o {output} -t 24
        """

rule vgcall:
    input:
        ref = "/data/wangj/myan/grape/17.callSVusingUpdateRef/VHP.hap1.fna",
        r1 = "/data/wangj/myan/grape/00.illumina.reads/{sample}_R1.fq.gz",
        r2 = "/data/wangj/myan/grape/00.illumina.reads/{sample}_R2.fq.gz",
        pack = rules.vgpack.output,
        gbz = rules.vgautoindex.output.gbz,
        min = rules.vgautoindex.output.min,
        dist = rules.vgautoindex.output.dist
    output:
        "08.vcf/{sample}.vcf"
    threads:
        32
    resources:
        mem_mb = 64000
    params:
        "{sample}"
    shell:
        """
        /data/wangj/myan/biotools/vg call \
        {input.gbz} -k {input.pack} -a -t 24 > {output}
        """
