SAMPLES, = glob_wildcards(config["input_folder"] + "/" + "{sample}.fa")

rule all:
    input:
        expand(config["output_folder"] + "/" + "{sample}.nucdiv",sample=SAMPLES)


rule axt2fa:
    input:
        config["input_folder"] + "/" + "{sample}.fa"
    output:
        config["output_folder"] + "/" + "{sample}.nucdiv"
    threads:
        2
    shell:
        """
        Rscript calculateNucDiv.R {input} {output}
        """