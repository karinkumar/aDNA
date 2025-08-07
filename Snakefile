configfile: "config.yaml"

rule all:
     input: expand("{sample}.sam", sample = config["samples"])


rule align:
     input: "{sample}.fastq.gz"
     output: "{sample}.sam"
     shell: """
~/software/bwa-mem2/bwa-mem2 mem -t 56 ~/hg38.fasta {input} > {output}
"""