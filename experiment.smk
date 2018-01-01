EXPR_PARAM, =  glob_wildcards("matrix/{expr_param}.csv")
EXPR_PARAM = EXPR_PARAM[:2]


METHODS_GENE = ['CC2Stat', 'DESeq2', 'sgRSEA', 'PBNPA']
METHODS_sgRNA = ['CC2Stat', 'DESeq2', 'edgeR']

rule all:
	input:
		expand("results/{method}/{expr_param}/gene.csv", method=METHODS_GENE, expr_param=EXPR_PARAM) +
		expand("results/{method}/{expr_param}/sgRNA.csv", method=METHODS_sgRNA, expr_param=EXPR_PARAM) +
		expand("results/mageck/{expr_param}.gene_summary.txt", expr_param=EXPR_PARAM) +
		expand("results/mageck/{expr_param}.sgrna_summary.txt", expr_param=EXPR_PARAM) 


rule run_mageck:
	input:
		"matrix/{expr_param}.mageck"
	params:
		c = "L1,L2,L3,L4",
		t = "H1,H2,H3,H4",
		dir_path = "results/mageck/{expr_param}"
	output:
		gene = "results/mageck/{expr_param}.gene_summary.txt",
		sgRNA = "results/mageck/{expr_param}.sgrna_summary.txt"
	benchmark:
		"benchmarks/mageck/{expr_param}.txt"
	shell:
		"mageck test -k {input} -t {params.t} -c {params.c} -n {params.dir_path}"
	

rule run_CC2Stat:
	input:
		"matrix/{expr_param}.csv"
	params:
		c = "L1,L2,L3,L4",
		t = "H1,H2,H3,H4"
	output:
		dir_path = "results/CC2Stat/{expr_param}",
		gene = "results/CC2Stat/{expr_param}/gene.csv",
		sgRNA = "results/CC2Stat/{expr_param}/sgRNA.csv"
	benchmark:
		"benchmarks/CC2Stat/{expr_param}.txt"
	shell:
		"CC2Stat.py -i {input} -o {output.dir_path} -c {params.t} -t {params.c}"


rule run_DESeq2:
	input:
		"matrix/{expr_param}.csv"
	params:
		m = "DESeq2"
	output:
		dir_path = "results/DESeq2/{expr_param}",
		gene = "results/DESeq2/{expr_param}/gene.csv",
		sgRNA = "results/DESeq2/{expr_param}/sgRNA.csv"
	benchmark:
		"benchmarks/DESeq2/{expr_param}.txt"
	shell:
		"Rscript --vanilla run_rtest.R {params.m} {input} {output.dir_path}"
		
	
rule run_sgRSEA:
	input:
		"matrix/{expr_param}.csv"
	params:
		m = "sgRSEA"
	output:
		dir_path = "results/sgRSEA/{expr_param}",
		gene = "results/sgRSEA/{expr_param}/gene.csv"
	benchmark:
		"benchmarks/sgRSEA/{expr_param}.txt"
	shell:
		"Rscript --vanilla run_rtest.R {params.m} {input} {output.dir_path}"

		
rule run_PBNPA:
	input:
		"matrix/{expr_param}.csv"
	params:
		m = "PBNPA"
	output:
		dir_path = "results/PBNPA/{expr_param}",
		gene = "results/PBNPA/{expr_param}/gene.csv"
	benchmark:
		"benchmarks/PBNPA/{expr_param}.txt"
	shell:
		"Rscript --vanilla run_rtest.R {params.m} {input} {output.dir_path}"

		
	
rule run_edgeR:
	input:
		"matrix/{expr_param}.csv"
	params:
		m = "edgeR"
	output:
		dir_path = "results/edgeR/{expr_param}",
		sgRNA = "results/edgeR/{expr_param}/sgRNA.csv"
	benchmark:
		"benchmarks/edgeR/{expr_param}.txt"
	shell:
		"Rscript --vanilla run_rtest.R {params.m} {input} {output.dir_path}"

		
	
