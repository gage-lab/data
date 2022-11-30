from snakemake.remote import AUTO
from snakemake.remote import FTP, HTTP

FTP = FTP.RemoteProvider()
HTTP = HTTP.RemoteProvider()


rule txome_gtf:
    input:
        FTP.remote(
            "ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.chr_patch_hapl_scaff.basic.annotation.gtf.gz",
            static=True,
        ),
    output:
        "ref/txome.chr{chrom}.gtf",
    shell:
        "zgrep -P ^chr{wildcards.chrom} {input} > {output}"


rule rmsk_gtf:
    input:
        HTTP.remote(
            "https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/GRCh38_GENCODE_rmsk_TE.gtf.gz",
            static=True,
        ),
    output:
        "ref/rmsk.chr{chrom}.gtf",
    shell:
        "zgrep -P ^chr{wildcards.chrom} {input} > {output}"


rule genome_fa:
    input:
        FTP.remote(
            "ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz",
            static=True,
        ),
    output:
        "ref/genome.chr{chrom}.fa",
    shell:
        "gzip -dc {input} > {output}"


rule txome_fa:
    input:
        fa=FTP.remote(
            "ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.transcripts.fa.gz",
            static=True,
        ),
        gtf=rules.txome_gtf.output,
    output:
        "ref/txome.chr{chrom}.fa",
    conda:
        "envs/seqkit.yaml"
    shell:
        """
        grep -o 'ENST[0-9]*\.[0-9]' {input.gtf} | sort | uniq | awk '{{print $1".*"}}' > names.lst
        zcat {input.fa} | seqkit grep -f names.lst -r -I > {output}
        rm -f names.lst
        """


rule reads:
    input:
        FTP.remote(
            [
                "ftp.ebi.ac.uk/biostudies/fire/E-GEUV-/001/E-GEUV-1/Files/E-GEUV-1/processed/NA20778.4.M_120208_1.bam",
                "ftp.ebi.ac.uk/biostudies/fire/E-GEUV-/001/E-GEUV-1/Files/E-GEUV-1/processed/NA20778.4.M_120208_1.bam.bai",
            ],
            static=True,
            keep_local=True,
        ),
    output:
        bam="bam/{sample}.chr{chrom}.bam",
        fq1="reads/{sample}.chr{chrom}.1.fq.gz",
        fq2="reads/{sample}.chr{chrom}.2.fq.gz",
    params:
        seed=lambda wildcards: abs(hash(wildcards.sample)) % 10000,
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools view -b -s{params.seed}.2 {input[0]} chr{wildcards.chrom} > {output.bam}
        samtools fastq -1 {output.fq1} -2 {output.fq2} {output.bam}
        """


rule all:
    input:
        expand(
            [
                "ref/txome.chr{chrom}.gtf",
                "ref/genome.chr{chrom}.fa",
                "ref/txome.chr{chrom}.fa",
                "ref/rmsk.chr{chrom}.gtf",
            ],
            chrom=21,
        ),
        expand(
            rules.reads.output,
            sample=["a", "b"],
            chrom=21,
        ),
