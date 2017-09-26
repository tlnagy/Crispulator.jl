all:	generate_file.jl convert_to_datamatrix.R
	julia generate_file.jl
	Rscript convert_to_datamatrix.R
	git add matrix/*
