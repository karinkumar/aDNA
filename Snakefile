configfile: "config.yaml"

rule all:
     input: expand("{sample}.sorted.bam", sample = config["samples"])


rule align:
     input: "{sample}.fastq.gz"
     output: "{sample}.sam"
     shell: """
~/software/bwa-mem2/bwa-mem2 mem -t 56 ~/hg38.fasta {input} > {output}
"""


rule bam:
     input: "{sample}.sam"
     output: "{sample}.sorted.bam"
     shell: """
     samtools view -bS {input} | samtools sort -o {output}
     samtools coverage {output} > {wildcards.sample}.cov.txt
     samtools index {output}
"""