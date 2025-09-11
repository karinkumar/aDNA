configfile: "config.yaml"

rule all:
     input: expand("{sample}_chr20_{cov}x.bam", sample = config["samples"], cov = config["COV"]), expand("gl/bcftoolsgenogvcfs{cov}x.vcf.gz", cov = config ["COV"]), expand("../val/{sample}_validation.vcf.gz", sample = config["samples"])


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


rule val_gt:
     input: "{sample}.sorted.bam"
     output: "../val/{sample}_validation.vcf.gz"
     params:
        sites="~/1kg30xASW/ref/1000GP.chr20.clean.nosingletons.sites.vcf.gz",
	sites_tsv="~/1kg30xASW/ref/1000GP.chr20.clean.nosingletons.sites.tsv.gz",
	mask="~/20160622.allChr.mask.bed",
	repeats="~/repeatsmask.bed",
	m_header="~/mask_header.txt",
	r_header="~/repeats_header.txt"
     shell:"""
OUT=/net/fantasia/home/kiranhk/aDNA/data/val/{wildcards.sample}.raw.Q20.q30.calls.vcf.gz
bcftools mpileup -f ~/hg38.fasta -I -E -a 'FORMAT/DP' --ignore-RG -T {params.sites} -Q 20 -q 30 -C 50 -r chr20 {input} | bcftools call -Aim -C alleles -T {params.sites_tsv} -Oz -o $OUT
#bcftools reheader -S {wildcards.sample} -o {output}
#remove masked sites
bcftools annotate -a {params.mask} -m +MASK=strict -h {params.m_header} -c CHROM,POS,FROM,TO $OUT | bcftools view -i 'INFO/MASK="strict"' | bcftools annotate -a {params.repeats} -m REPEATS=repeats -h {params.r_header} -c CHROM,POS,FROM,TO | bcftools view --exclude 'REPEATS="repeats"' -Oz -o ../val/{wildcards.sample}_removedmaskrepeats.vcf.gz
bcftools index -f  ../val/{wildcards.sample}_removedmaskrepeats.vcf.gz

#depth and QUAL filter
DOC=$(bcftools query -f '%INFO/DP\n' $OUT |  awk 'BEGIN {{ s = 0; l=0; }} {{ s+=$1; l++; }} END {{ print s/l;}}')

LOW=$(python3 limitsDoC.py $DOC | awk '{{print $1}}')
UPP=$(python3 limitsDoC.py $DOC | awk '{{print $2}}')
echo $LOW $UPP

bcftools filter --exclude "FORMAT/DP<$LOW |  FMT/DP>$UPP & QUAL<30" ../val/{wildcards.sample}_removedmaskrepeats.vcf.gz -Oz -o {output}
bcftools index -f {output}
"""


#rule combine_valgt:
     #input: "../val/{sample}_validation.vcf.gz"
     #output: "../val/cleaned_validation.vcf.gz"
     #shell:"""
#combine files


