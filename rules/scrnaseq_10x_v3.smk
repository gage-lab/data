# download reference data from 10x Genomics
rule get_refdata:
    input:
        HTTP.remote(
            "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz",
            static=True,
        ),
    output:
        multiext("scrnaseq_10x_v3/ref/", "genome.chr{chrom}.fa", "genes.chr{chrom}.gtf"),
    cache: "omit-software"
    shell:
        """
        tar -xf {input} --wildcards '*genes.gtf' '*genome.fa'
        seqkit grep -p "chr{wildcards.chrom}" $(basename {input} .tar.gz)/fasta/genome.fa > {output[0]}
        grep "^chr{wildcards.chrom}" $(basename {input} .tar.gz)/genes/genes.gtf > {output[1]}
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
        "scrnaseq_10x_v3/10x_v3_whitelist.txt",
    cache: "omit-software"
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
        "scrnaseq_10x_v3/{sample}.chr{chrom}.ids.txt",
    params:
        seed=lambda wildcards: abs(hash(wildcards.sample)) % 10000,
        subsample=0.5,
    conda:
        "../environment.yaml"
    shell:
        """
        touch -m {input[1]}
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
    cache: "omit-software"
    conda:
        "../environment.yaml"
    shell:
        "tar -xf {input.fastq}"


rule get_reads:
    input:
        untar=rules.untar_reads.output,
        ids=rules.find_rids.output,
    output:
        "scrnaseq_10x_v3/{sample}.chr{chrom}_{read}.fastq.gz",
    conda:
        "../environment.yaml"
    shell:
        "seqkit grep -n -f {input.ids} {input.untar}/*{wildcards.read}_001.fastq.gz > {output}"


rule scrnaseq_10x_v3:
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
