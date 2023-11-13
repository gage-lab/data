# ONT from https://github.com/nanopore-wgs-consortium/NA12878/blob/master/RNA.md
# chose those files failing QC as tehy are smaller
def get_reads(wildcards):
    if wildcards.libtype == "ONT_directRNA":
        return HTTP.remote(
            [
                "https://s3.amazonaws.com/nanopore-human-wgs/rna/bamFiles/NA12878-DirectRNA.fail.dedup.NoU.fastq.hg38.minimap2.sorted.bam",
                "https://s3.amazonaws.com/nanopore-human-wgs/rna/bamFiles/NA12878-DirectRNA.fail.dedup.NoU.fastq.hg38.minimap2.sorted.bam.bai",
            ],
            static=True,
            keep_local=True,
        )
    elif wildcards.libtype == "ONT_cDNA":
        return HTTP.remote(
            [
                "https://s3.amazonaws.com/nanopore-human-wgs/rna/bamFiles/NA12878-cDNA-1D.fail.dedup.fastq.hg38.minimap2.sorted.bam",
                "https://s3.amazonaws.com/nanopore-human-wgs/rna/bamFiles/NA12878-cDNA-1D.fail.dedup.fastq.hg38.minimap2.sorted.bam.bai",
            ],
            static=True,
            keep_local=True,
        )
    else:
        raise ValueError("Unknown libtype: {}".format(wildcards.libtype))


rule longrnaseq_reads:
    input:
        get_reads,
    output:
        bam="longrnaseq/{libtype}/{sample}.chr{chrom}.bam",
        fq="longrnaseq/{libtype}/{sample}.chr{chrom}.fq.gz",
    conda:
        "../environment.yaml"
    shell:
        """
        samtools view -b {input[0]} chr{wildcards.chrom} > {output.bam}
        samtools fastq {output.bam} | gzip > {output.fq}
        """


rule longrnaseq:
    input:
        expand(
            rules.longrnaseq_reads.output,
            sample=["a", "b"],
            libtype=["ONT_directRNA", "ONT_cDNA"],
            chrom=21,
        ),
