include(joinpath(Pkg.dir("Crispulator"), "src", "simulation", "load.jl"))


srand(577322681)

num_sample = 4
for seq_depth in [ 10, 100, 1000 ]
    for bin_prob in [0.1, 0.25]
        for σ in [0.01, 0.1, 0.2]
            for pheno_prob in [0, 0.01, 0.2]
                facs_param = FacsScreen()
                facs_param.num_genes = 1000
                facs_param.coverage = 10
                facs_param.representation = 100
                facs_param.bin_info[:bin1] = (0, bin_prob)
                facs_param.bin_info[:bin2] = (1-bin_prob, 1)
                facs_param.seq_depth = seq_depth
                facs_param.σ = σ


                lib = Library(CRISPRn())
                lib.phenotype_probs = Categorical([pheno_prob/2,0.05,0.95-pheno_prob,pheno_prob/2])
                guides, guide_freqs_dist = construct_library(facs_param, lib)

                dir_path = string("data/facs_",  seq_depth, "_", @sprintf("%.2f", bin_prob), "_", @sprintf("%.2f",σ), "_", @sprintf("%.2f",pheno_prob), "/")
                println(dir_path)
                mkpath(dir_path)
                for rep_no in 1:num_sample
                    println("Running $(rep_no)...")
                    cells, cell_phenotypes = transfect(facs_param, lib, guides, guide_freqs_dist)
                    bin_cells = select(facs_param, cells, cell_phenotypes, guides)
                    freqs = counts_to_freqs(bin_cells, length(guides))

                    seq_depths = Dict{Symbol, Int}()
                    for binname in keys(bin_cells)
                        seq_depths[binname] = rand(Poisson(facs_param.seq_depth))
                    end

                    raw_data = sequencing(seq_depths, guides, freqs)
                    bc_counts, genes = differences_between_bins(raw_data)

                    file_name = string(dir_path,"/bc_count_$(rep_no).csv")
                    writetable(file_name, bc_counts)
                end
            end
        end
    end
end
