SAMPLES, = glob_wildcards(config["input_folder"] + "/" + "{sample}.axt")

rule all:
    input:
        expand(config["output_folder"] + "/" + "{sample}.kaks",sample=SAMPLES)


rule axt2fa:
    input:
        config["input_folder"] + "/" + "{sample}.axt"
    output:
        config["output_folder"] + "/" + "{sample}.kaks"
    threads:
        2
    shell:
        """
        /home/myan/biotools/software.package/KaKs_Calculator3.0/src/KaKs -i \
        {input} -o {output} -c 1 -m YN
        """