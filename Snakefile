storage ftp:
	provider="ftp"
storage http:
	provider="http"

conda:
	"envs/global.yaml"

include: "rules/rnaseq.smk"
include: "rules/scrnaseq_10x_v3.smk"
include: "rules/wgs.smk"
include: "rules/longrnaseq.smk"
