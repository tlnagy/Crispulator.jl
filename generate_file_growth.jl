include(joinpath(Pkg.dir("Crispulator"), "src", "simulation", "load.jl"))

srand(577322681)

num_sample = 4
for seq_depth in [ 10, 100, 1000 ]
    for num_bottlenecks in [1, 10]
            for pheno_prob in [0, 0.01, 0.2]
                for err in [0.01, 0.1, 0.2]
                    growth_param = GrowthScreen()
                    growth_param.num_genes = 1000
                    growth_param.coverage = 10
                    growth_param.representation = 100
                    growth_param.seq_depth = seq_depth
                    growth_param.noise = err
                    max_phenotype_dists = Dict{Symbol, Tuple{Float64, Sampleable}}(
                        :inactive => (1-0.05-pheno_prob, Delta(0.0)),
                        :negcontrol => (0.05, Delta(0.0)),
                        :increasing => (pheno_prob/2, TruncatedNormal(0.1, 0.1, 0.025, 1)),
                        :decreasing => (pheno_prob/2, TruncatedNormal(-0.55, 0.2, -1, -0.1))
                    )

                    lib = Library(max_phenotype_dists, CRISPRn())
                    guides, guide_freqs_dist = construct_library(growth_param, lib)

                    dir_path = string("data/growth_",  seq_depth, "_", @sprintf("%d", num_bottlenecks), "_", @sprintf("%.2f",err), "_", @sprintf("%.2f",pheno_prob), "/")
                    println(dir_path)
                    mkpath(dir_path)
                    for rep_no in 1:num_sample
                        println("Running $(rep_no)...")
                        cells, cell_phenotypes = transfect(growth_param, lib, guides, guide_freqs_dist)
                        bin_cells = select(growth_param, cells, cell_phenotypes, guides)
                        freqs = counts_to_freqs(bin_cells, length(guides))

                        seq_depths = Dict{Symbol, Int}()
                        for binname in keys(bin_cells)
                            seq_depths[binname] = rand(Poisson(growth_param.seq_depth))
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
