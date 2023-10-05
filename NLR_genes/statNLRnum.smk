SAMPLES, = glob_wildcards("03.candidate.NLR.peps.Ortho/{sample}_nlr.fa")

print("Total sample size: ",len(SAMPLES))

rule all:
    input:
        expand("04.NLR.num/{sample}.txt",sample=SAMPLES)


rule statNLRnum:
    input:
        "03.candidate.NLR.peps.Ortho/{sample}_nlr.fa"
    output:
        "04.NLR.num/{sample}.txt"
    threads:
        4
    resources:
        mem_mb = 4000
    params:
        "{sample}"
    run:
        from Bio import SeqIO
        fw = open(output[0],'w')
        i = 0
        for rec in SeqIO.parse(input[0],'fasta'):
            i += 1
        fw.write("%s\t%d\n"%(params[0],i))