# download genome fasta
rule wgs_genome_fa:
    input:
        FTP.remote(
            "ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz",
            static=True,
            keep_local=True,
            immediate_close=True,
            timeout=600,
        ),
    output:
        "wgs/ref/genome.chr20.fa",
    cache: "omit-software"
    conda:
        "../environment.yaml"
    shell:
        "gzip -dc {input} | seqkit grep -n -p chr20 > {output}"


rule wgs_reads:
    input:
        FTP.remote(
            [
                "ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/phase3/data/NA20778/alignment/NA20778.chrom20.ILLUMINA.bwa.TSI.low_coverage.20130415.bam",
                "ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/phase3/data/NA20778/alignment/NA20778.chrom20.ILLUMINA.bwa.TSI.low_coverage.20130415.bam.bai",
            ],
            static=True,
            keep_local=True,
            immediate_close=True,
            timeout=600,
        ),
    output:
        bam="wgs/{sample}.chr20.bam",
        fq1="wgs/{sample}.chr20.1.fq.gz",
        fq2="wgs/{sample}.chr20.2.fq.gz",
    params:
        seed=lambda wildcards: abs(hash(wildcards.sample)) % 10000,
        subsample=0.4,
    conda:
        "../environment.yaml"
    shell:
        """
        samtools view -b --subsample-seed {params.seed} --subsample {params.subsample} {input} > {output.bam}
        samtools fastq -1 {output.fq1} -2 {output.fq2} -0 /dev/null -s /dev/null {output.bam}
        """


rule wgs:
    input:
        expand(
            rules.wgs_genome_fa.output,
        ),
        expand(
            rules.wgs_reads.output,
            sample=["a"],
        ),
