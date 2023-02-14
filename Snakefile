from snakemake.remote import FTP, HTTP

FTP = FTP.RemoteProvider()
HTTP = HTTP.RemoteProvider()


include: "rules/rnaseq.smk"
include: "rules/scrnaseq_10x_v3.smk"
include: "rules/wgs.smk"
