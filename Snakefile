from snakemake.remote import FTP, HTTP

FTP = FTP.RemoteProvider()
HTTP = HTTP.RemoteProvider()


include: "rules/rnaseq.smk"
include: "rules/dsc_rnaseq.smk"
include: "rules/wgs.smk"
