# download reference data from 10x Genomics
rule untar_refdata:
    input:
        HTTP.remote(
            "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz",
            static=True,
        ),
    output:
        directory("refdata-gex-GRCh38-2020-A"),
    shell:
        "tar -xf {input}"


rule get_refdata:
    input:
        rules.untar_refdata.output,
    output:
        fa="dsc_rnaseq/ref/genome.chr{chrom}.fa",
        gtf="dsc_rnaseq/ref/genes.chr{chrom}.gtf",
    conda:
        "../environment.yaml"
    shell:
        """
        seqkit grep -n -p "chr{wildcards.chrom}" {input}/fasta/genome.fa > {output.fa}
        grep "^chr{wildcards.chrom}" {input}/genes/genes.gtf > {output.gtf}
        """


# download barcode whitelist from 10x Genomics
# for 10x v3 chemistry
rule get_whitelist:
    input:
        HTTP.remote(
            "https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz",
            static=True,
        ),
    output:
        "dsc_rnaseq/10x_v3_whitelist.txt",
    shell:
        "gunzip -c {input} > {output}"


rule find_rids:
    input:
        HTTP.remote(
            [
                "https://cf.10xgenomics.com/samples/cell-exp/6.0.0/Brain_Tumor_3p_LT/Brain_Tumor_3p_LT_possorted_genome_bam.bam",
                "https://cf.10xgenomics.com/samples/cell-exp/6.0.0/Brain_Tumor_3p_LT/Brain_Tumor_3p_LT_possorted_genome_bam.bam.bai",
            ],
            static=True,
            keep_local=True,
        ),
    output:
        "dsc_rnaseq/{sample}.chr{chrom}.ids.txt",
    params:
        seed=lambda wildcards: abs(hash(wildcards.sample)) % 10000,
        subsample=0.5,
    conda:
        "../environment.yaml"
    shell:
        """
        # get read IDs from BAM file
        samtools view {input[0]} chr{wildcards.chrom} | \
            awk 'BEGIN {{srand({params.seed})}} {{if (rand() < {params.subsample}) print $1}}' > {output}
        """


rule untar_reads:
    input:
        fastq=HTTP.remote(
            "https://cf.10xgenomics.com/samples/cell-exp/6.0.0/Brain_Tumor_3p_LT/Brain_Tumor_3p_LT_fastqs.tar",
            static=True,
        ),
    output:
        directory("Brain_Tumor_3p_LT_fastqs"),
    conda:
        "../environment.yaml"
    shell:
        "tar -xf {input.fastq}"


rule get_reads:
    input:
        untar=rules.untar_reads.output,
        ids=rules.find_rids.output,
    output:
        "dsc_rnaseq/{sample}.chr{chrom}_{read}.fastq.gz",
    conda:
        "../environment.yaml"
    shell:
        "seqkit grep -n -f {input.ids} {input.untar}/*{wildcards.read}_001.fastq.gz > {output}"


rule dsc_rnaseq:
    input:
        expand(
            rules.get_reads.output,
            sample=["a", "b"],
            chrom=21,
            read=["I1", "I2", "R1", "R2"],
        ),
        expand(
            rules.get_refdata.output,
            chrom=21,
        ),
