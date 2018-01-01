NCORE = 4

all:	experiment.smk
	snakemake -s experiment.smk --cores $(NCORE)


force:	experiment.smk
	snakemake -s experiment.smk --forceall --cores $(NCORE)


dry:	experiment.smk
	snakemake -s experiment.smk -np
