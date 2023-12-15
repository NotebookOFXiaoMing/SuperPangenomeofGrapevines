SAMPLES, = glob_wildcards(config["input_folder"] + "/" + "{sample}.axt")

rule all:
    input:
        expand(config["output_folder"] + "/" + "{sample}.fa",sample=SAMPLES)


rule axt2fa:
    input:
        config["input_folder"] + "/" + "{sample}.axt"
    output:
        config["output_folder"] + "/" + "{sample}.fa"
    threads:
        2
    shell:
        """
        python axt2fasta.py {input} {output}
        """