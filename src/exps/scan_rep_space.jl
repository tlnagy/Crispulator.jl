# There are many 3 primary stages of the experiment where representation
# is important: at transfection, at the bottleneck(s), and at
# sequencing.

function scan_rep_space(filepath; debug=false, quiet=false)

    if !debug
        parameters = Dict{Symbol, Vector}(
            :representation => map(x->round(Int, x), logspace(0, 4, 13)),
            :bottleneck_representation => map(x->round(Int, x),  logspace(0,4,13)),
            :seq_depth => map(x->round(Int, x), logspace(0, 4, 13))
        )
        num_runs = 10
    else
        parameters = Dict{Symbol, Vector}(
            :representation => map(x->round(Int, x), logspace(0, 2, 2)),
            :bottleneck_representation => map(x->round(Int, x),  logspace(0,2,2)),
            :seq_depth => [10^2]
        )
        num_runs = 1
    end

    const overlap = intersect(fieldnames(FacsScreen), fieldnames(GrowthScreen))
    # custom function for handling both growth and a facs screen in the
    # same relational datastructure
    flatten_overlap = (setup, lib) -> begin
        local results = Any[typeof(setup)]
        for name in overlap
            push!(results, getfield(setup, name))
        end
        results
    end

    methods = [venn, auprc]
    measures = [:inc, :dec, :incdec]
    genetypes = [:sigmoidal, :linear, :all]
    test_method_wrapper = (bc_counts, genes) -> test_methods(genes, methods, measures, genetypes)

    results = []
    for screentype in [FacsScreen(), GrowthScreen()]
        for crisprtype in [CRISPRi(), CRISPRn()]
            if typeof(screentype) == GrowthScreen
                max_phenotype_dists = Dict{Symbol, Tuple{Float64, Sampleable}}(
                    :inactive => (0.83, Delta(0.0)),
                    :negcontrol => (0.05, Delta(0.0)),
                    :increasing => (0.02, TruncatedNormal(0.1, 0.1, 0.025, 1)),
                    :decreasing => (0.1, TruncatedNormal(-0.55, 0.2, -1, -0.1))
                )
                lib = Library(max_phenotype_dists, crisprtype)
            else
                lib = Library(crisprtype)
            end

            runs = build_parameter_space(screentype, parameters, num_runs)

            before = time()
            result = pmap(args -> run_exp(args[1], lib, test_method_wrapper; run_idx=args[2], flatten_func=flatten_overlap), runs)
            (!quiet) && println("$(time() - before) seconds")
            result = DataFrame(permutedims(hcat(result...), [2, 1]))
            result[:crisprtype] = typeof(crisprtype)
            push!(results, result)
        end
    end
    results = vcat(results...)
    hierarchy = vcat([hcat(item...) for item in IterTools.product(map(Symbol, methods), measures, genetypes)]...)
    new_names = [[:method, :measure, :genetype, :score, :screen]...; overlap...; :run; :crisprtype]
    results = construct_hierarchical_label(hierarchy, results, new_names)
    CSV.write(filepath, results)
end
