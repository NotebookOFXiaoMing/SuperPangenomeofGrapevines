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

