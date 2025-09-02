configfile: "config.yaml"

rule all:
     input: expand("{sample}_chr20_{cov}x.bam", sample = config["samples"], cov = config["COV"]), expand("gl/bcftoolsgenogvcfs{cov}x.vcf.gz", cov = config ["COV"])


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

rule downsample:
    input: "{sample}.sorted.bam"
    output:
        "{sample}_chr20_{cov}x.bam"
        #bai=expand("{sample}_chr20_1x.bam.bai", sample = sample_list)
    params:
        hg38 = "~/hg38.fasta",  # reference genome
        region = "chr20"
    shell:
        """
	samtools index {input}
        set -e  # stop on error
        COV=$(samtools coverage -r chr20 {input} | awk 'NR==2 {{print $7}}')
        echo $COV
        # Variables
        TARG_COV={wildcards.cov}
        FRAC=$(echo "scale=10; $TARG_COV/$COV" | bc)  # calculating fraction for downsampling, scale=10 ensures proper float calculation

        # Check the input file integrity
        samtools quickcheck {input}

        # Perform downsampling if output doesn't exist
        if [ ! -f {output} ]; then
            echo "Performing downsampling"
            samtools view -T {params.hg38} -s $FRAC -bo {output} {input} {params.region}
            samtools index {output}
        else
            echo "{output} already exists, skipping downsampling."
        fi

        # Check coverage
        """


#for gl insert correct gl sites that need to be updated anyway

rule gl:
    input:
        "{sample}_chr20_{cov}x.bam"
    output:
        vcf="{sample}_chr20.{cov}x.vcf.gz",
        tbi="{sample}_chr20.{cov}x.vcf.gz.csi"
    params:
        sites="~/1kg30xASW/ref/1000GP.chr20.clean.nosingletons.sites.vcf.gz",
	sites_tsv="~/1kg30xASW/ref/1000GP.chr20.clean.nosingletons.sites.tsv.gz"
    shell:
        """
        samtools index {input} #force reindexing to stay current
        ID=$(echo {input} | grep -o ^'[^.]*[.]')
        bcftools mpileup -f ~/hg38.fasta -I -E -a 'FORMAT/DP' -T {params.sites}  -r chr20 {input}  -Ou | bcftools call -Aim -C alleles -T {params.sites_tsv} -Oz -o {output.vcf}
        bcftools index {output.vcf}
        """


rule combine_gl:
    input: expand("{sample}_chr20.{cov}x.vcf.gz", sample = config["samples"], cov = config["COV"])
    output:
        "gl/bcftoolsgenogvcfs{cov}x.vcf.gz"
    shell:
        """
        ls ~/aDNA/data/*chr20.{wildcards.cov}x.vcf.gz > {wildcards.cov}list.txt
        bcftools merge -m none -r chr20 -Oz -o {output} -l {wildcards.cov}list.txt
        bcftools index {output}
        """
