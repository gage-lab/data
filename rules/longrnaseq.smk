rule long_rna_reads:
    output:
        bam="longrnaseq/{sample}.chr{chrom}.bam",
        fq="longrnaseq/{sample}.chr{chrom}.fq.gz",
    conda:
        "../environment.yaml"
    shell:
        """
        aws s3 --no-sign-request sync s3://sg-nex-data/data/sequencing_data_ont/bam/genome/SGNex_MCF7_directcDNA_replicate1_run2/SGNex_MCF7_directcDNA_replicate1_run2.bam .
        aws s3 --no-sign-request sync s3://sg-nex-data/data/sequencing_data_ont/bam/genome/SGNex_MCF7_directcDNA_replicate1_run2/SGNex_MCF7_directcDNA_replicate1_run2.bam.bai .
        samtools view -b SGNex_MCF7_directcDNA_replicate1_run2.bam {wildcards.chrom} > {output.bam}
        samtools fastq {output.bam} | gzip > {output.fq}
        """


rule longrnaseq:
    input:
        expand(
            rules.reads.output,
            sample=["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
            chrom=21,
        ),
