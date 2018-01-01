all:	generate_file_facs.jl generate_file_growth.jl convert_to_datamatrix.R
	rm -rf matrix/*
	rm -rf data/*	
	julia generate_file_facs.jl
	julia generate_file_growth.jl
	Rscript convert_to_datamatrix.R
	git add matrix/*
