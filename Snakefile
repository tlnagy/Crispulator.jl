seq_depth = [ 10, 100, 1000 ]
rep = [ 10, 100, 1000 ]
pheno = [ 0.01, 0.1, 0.2 ]
err = [ 0.01, 0.1, 0.2 ]
bin_ratio = [ 0.1, 0.25 ]
random_seed = 10005
from util import *

rule all:
	input:
		expand("data/FACS_{seed}_{br}_{r}_{s}_{p}_{e}", seed=random_seed, r=rep, s=seq_depth, p=pheno, e=err, br=bin_ratio) +
		expand("data/GROWTH_{seed}_{r}_{s}_{p}_{e}", seed=random_seed, r=rep, s=seq_depth, p=pheno, e=err) +
		expand("data/FACS_{seed}_{br}_{r}_{s}_{p}_{e}/bc_count_4.csv", seed=random_seed, r=rep, s=seq_depth, p=pheno, e=err, br=bin_ratio) +
		expand("data/GROWTH_{seed}_{r}_{s}_{p}_{e}/bc_count_4.csv", seed=random_seed, r=rep, s=seq_depth, p=pheno, e=err) +
		expand("matrix/FACS_{seed}_{br}_{r}_{s}_{p}_{e}.csv", seed=random_seed, r=rep, s=seq_depth, p=pheno, e=err, br=bin_ratio) +
		expand("matrix/GROWTH_{seed}_{r}_{s}_{p}_{e}.csv", seed=random_seed, r=rep, s=seq_depth, p=pheno, e=err) 


rule merge_facs:
	input:
		in_dir = "data/FACS_{seed}_{br}_{r}_{s}_{p}_{e}",
		in_file = "data/FACS_{seed}_{br}_{r}_{s}_{p}_{e}/bc_count_4.csv"
	output:
		out_file = "matrix/FACS_{seed}_{br}_{r}_{s}_{p}_{e}.csv"
	run:
		merge_bc_counts(str(input.in_dir), str(output.out_file))
	

rule merge_growth:
	input:
		in_dir = "data/GROWTH_{seed}_{r}_{s}_{p}_{e}",
		in_file = "data/GROWTH_{seed}_{r}_{s}_{p}_{e}/bc_count_4.csv"
	output:
		out_file = "matrix/GROWTH_{seed}_{r}_{s}_{p}_{e}.csv"
	run:
		merge_bc_counts(str(input.in_dir), str(output.out_file))
	
	
		
rule run_facs:
	input:
		"simdata_generator.jl"
	output:
		out_dir = "data/FACS_{seed}_{br}_{r}_{s}_{p}_{e}",
		out_file = "data/FACS_{seed}_{br}_{r}_{s}_{p}_{e}/bc_count_4.csv"

	shell:
		"julia {input} --screen_type FACS --representation {wildcards.r} " +
		"--seq_depth {wildcards.s} " +
		"--positive_ratio {wildcards.p} " +
		"--noise_ratio {wildcards.e} " +
		"--bin_ratio {wildcards.br} " +
		"--output_path {output.out_dir}"


rule run_growth:
	input:
		"simdata_generator.jl"
	output:
		out_dir = "data/GROWTH_{seed}_{r}_{s}_{p}_{e}",
		out_file = "data/GROWTH_{seed}_{r}_{s}_{p}_{e}/bc_count_4.csv"

	shell:
		"julia {input} --screen_type GROWTH --representation {wildcards.r} " +
		"--seq_depth {wildcards.s} " +
		"--positive_ratio {wildcards.p} " +
		"--noise_ratio {wildcards.e} " +
		"--output_path {output.out_dir}"

