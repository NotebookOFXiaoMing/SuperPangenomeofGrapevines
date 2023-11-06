SAMPLES, = glob_wildcards(config["qrys"])

print(SAMPLES)


rule all:
    input:
        expand("02.syri.output/{sample}/{sample}_syri.vcf",sample=SAMPLES),
        expand("02.syri.output/{sample}/{sample}.pdf",sample=SAMPLES),
        expand("03.InvTrans.vcf/{sample}_InvTrans.vcf",sample=SAMPLES),
        expand("03.InDelLarger50.vcf/{sample}_InsDel.vcf",sample=SAMPLES),
        expand("04.InvTrans.len/{sample}.txt",sample=SAMPLES),
        expand("04.InDel.len/{sample}.txt",sample=SAMPLES),
        "05.SvLenandNum/SvLenandNum.txt",
        "05.SvLenandNum/SvLenandNum.pdf",
        "06.merged.vcf/merged.InsDel.vcf",
        "06.merged.vcf/merged.InsDel.filter.vcf",
        "06.merged.vcf/merged.InsDel.filter.sorted.vcf.gz",
        config['graph_prefix'] + "giraffe.gbz",


rule minimap2:
    input:
        ref = config['ref'],
        qry = config["qrys"]
    output:
        "01.sam/{sample}.sam"
    threads:
        16
    resources:
        mem_mb = 48000
    shell:
        """
        minimap2 -ax asm5 --eqx {input.ref} {input.qry} -o {output} -t {threads}
        """

rule syri_call:
    input:
        ref = config['ref'],
        qry = config["qrys"],
        sam = rules.minimap2.output
    output:
        vcf = "02.syri.output/{sample}/{sample}_syri.vcf",
        txt = "02.syri.output/{sample}/{sample}_syri.out"
    threads:
        16
    resources:
        mem_mb = 24000
    params:
        output_folder = "02.syri.output/{sample}",
        prefix = "{sample}_"
    shell:
        """
        syri -c {input.sam} -r {input.ref} -q {input.qry} \
        -k -F S --nc {threads} --dir {params.output_folder} \
        --prefix {params.prefix}
        """


rule vim_genomes:
    input:
        ref = config['ref'],
        qry = config['qrys'],
        vcf = rules.syri_call.output
    output:
        "02.syri.output/{sample}/genomes.txt"
    threads:
        2
    resources:
        mem_mb = 8000
    params:
        ref = "ref",
        qry = "{sample}"
        
    run:
        fw = open(output[0],'w')
        fw.write("%s\t%s\n%s\t%s\n"%(input[0],params[0],input[1],params[1]))
        fw.close()

rule plotsr:
    input:
        syri_out = rules.syri_call.output.txt,
        genomes_txt = rules.vim_genomes.output
    output:
        "02.syri.output/{sample}/{sample}.pdf"
    threads:
        8
    resources:
        mem_mb = 16000
    shell:
        """
        plotsr --sr {input.syri_out} --genomes {input.genomes_txt} -W 6 -H 4 -o {output}
        """

rule filterIndelLarger50:
    input:
        rules.syri_call.output.vcf
    output:
        "03.InDelLarger50.vcf/{sample}_InsDel.vcf",
        "04.InDel.len/{sample}.txt"
    threads:
        4
    resources:
        mem_mb = 8000
    params:
        "{sample}"
    shell:
        """
        python scripts/filterSyriVCFInDel.py {params} {input} {output[0]} {output[1]}
        """

rule filterInvTrans:
    input:
        rules.syri_call.output.vcf
    output:
        "03.InvTrans.vcf/{sample}_InvTrans.vcf",
        "04.InvTrans.len/{sample}.txt"
    threads:
        4
    resources:
        mem_mb = 8000
    params:
        "{sample}"
    shell:
        """
        python scripts/filterSyriVCFInvTrans.py {params} {input} {output[0]} {output[1]}
        """

rule LenAndNumber:
    input:
        a = expand("04.InvTrans.len/{sample}.txt",sample=SAMPLES),
        b = expand("04.InDel.len/{sample}.txt",sample=SAMPLES)
    output:
        "05.SvLenandNum/SvLenandNum.txt"
    threads:
        2
    resources:
        mem_mb = 4000
    shell:
        """
        cat {input.a} {input.b} > {output}
        """

rule figures:
    input:
        rules.LenAndNumber.output
    output:
        "05.SvLenandNum/SvLenandNum.pdf"
    threads:
        4
    resources:
        mem_mb = 16000
    script:
        "scripts/figures.R"

rule mergeInsDel:
    input:
        expand("03.InDelLarger50.vcf/{sample}_InsDel.vcf",sample=SAMPLES)
    output:
        "06.merged.vcf/merged.InsDel.vcf"
    threads:
        16
    resources:
        mem_mb = 16000
    params:
        "03.InDelLarger50.vcf/*.vcf"
    shell:
        """
        ls {input} > vcf.list
        SURVIVOR merge vcf.list 1000 0 1 1 0 50 {output}
        """

rule filterMergedVCF:
    input:
        ref = config['ref'],
        fai = config['fai'],
        vcf = rules.mergeInsDel.output
    output:
        vcf = "06.merged.vcf/merged.InsDel.filter.vcf"
    threads:
        4
    resources:
        mem_mb = 8000
    shell:
        """
        python scripts/filterVCFandREFinconsistentLoci.py {input.vcf} {input.ref} {output.vcf} {input.fai}
        """


rule bcftools_sort:
    input:
        rules.filterMergedVCF.output.vcf
    output:
        vcf = "06.merged.vcf/merged.InsDel.filter.sorted.vcf.gz"
    threads:
        16
    resources:
        mem_mb = 24000
    shell:
        """
        /home/myan/biotools/software.package/bcftools-1.17/bcftools sort {input} -O z -o {output.vcf}
        tabix {output.vcf}
        """

rule vgautoindex:
    input:
        ref = config['ref'],
        vcf = rules.bcftools_sort.output.vcf
    output:
        gbz = config['graph_prefix'] + "giraffe.gbz",
        dist = config['graph_prefix'] + ".dist",
        min = config['graph_prefix'] + ".min"
    threads:
        24
    resources:
        mem_mb = 64000
    params:
        config['graph_prefix']
    shell:
        """
        mkdir tmpdir
        vg autoindex --workflow giraffe -r {input.ref} -v {input.vcf} -p {params} -t {threads} -T tmpdir/
        """
