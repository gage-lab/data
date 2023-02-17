# download genome fasta
rule wgs_genome_fa:
    input:
        FTP.remote(
            "ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa",
            static=True,
            keep_local=True,
            immediate_close=True,
        ),
    output:
        "wgs/ref/genome.chr{chrom}.fa",
    cache: "omit-software"
    conda:
        "../environment.yaml"
    shell:
        "seqkit grep -n -p chr{wildcards.chrom} {input} > {output}"


rule wgs_reads:
    output:
        cram="wgs/{sample}.chr{chrom}.cram",
        fq1="wgs/{sample}.chr{chrom}.1.fq.gz",
        fq2="wgs/{sample}.chr{chrom}.2.fq.gz",
    params:
        seed=lambda wildcards: abs(hash(wildcards.sample)) % 10000,
        subsample=0.4,
        cram="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/TSI/NA20778/alignment/NA20778.alt_bwamem_GRCh38DH.20150718.TSI.low_coverage.cram",
        ref="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa",
    conda:
        "../environment.yaml"
    shell:
        """
        samtools view --cram -T {params.ref} --subsample-seed {params.seed} --subsample {params.subsample} {params.cram} chr{wildcards.chrom} > {output.cram}
        samtools fastq -1 {output.fq1} -2 {output.fq2} -0 /dev/null -s /dev/null {output.cram}
        """


rule wgs:
    input:
        expand(
            rules.wgs_genome_fa.output,
            chrom=21,
        ),
        expand(
            rules.wgs_reads.output,
            sample=["a"],
            chrom=21,
        ),
