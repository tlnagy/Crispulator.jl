# This experiment will compare the behavior of pvalue, effect size, and
# product over a wide range of screen and library designs

function compare_methods(filepath; debug=false, quiet=false)

    runs = []

    if !debug
        num_runs = 25
        levels = (10, 100, 1000)
    else
        num_runs = 1
        levels = (10)
    end

    screen_combos = (((FacsScreen, NaN), (GrowthScreen, 10), (GrowthScreen, 20)),
                    levels,
                    (CRISPRi(), CRISPRn()),
                    (Linear(), Sigmoidal()))

    for (screentype, representation, crisprtype, genetype) in IterTools.product(screen_combos...)
        screen = screentype[1]()
        screen.representation = representation
        screen.bottleneck_representation = representation
        screen.seq_depth = representation
        if screentype[1] == FacsScreen
            screen.bin_info = OrderedDict{Symbol, Tuple{Float64, Float64}}(:bin1 => (0.0, 0.25), :bin2 => (0.75, 1.0))
            max_phenotype_dists = Dict{Symbol, Tuple{Float64, Sampleable}}(
                :inactive => (0.75, Delta(0.0)),
                :negcontrol => (0.05, Delta(0.0)),
                :increasing => (0.1, TruncatedNormal(0.55, 0.2, 0.1, 1)),
                :decreasing => (0.1, TruncatedNormal(-0.55, 0.2, -1, -0.1))
            )
        else
            screen.num_bottlenecks = screentype[2]
            max_phenotype_dists = Dict{Symbol, Tuple{Float64, Sampleable}}(
                :inactive => (0.85, Delta(0.0)),
                :negcontrol => (0.05, Delta(0.0)),
                :decreasing => (0.1, TruncatedNormal(-0.55, 0.2, -1, -0.1))
            )
        end

        kd_phenotype_relationships = Dict{Symbol, Tuple{Float64, KDPhenotypeRelationship}}(
            Symbol(typeof(genetype)) => (1.0, genetype)
        )

        if typeof(crisprtype) == CRISPRn
            knockdown_dist = Dict{Symbol, Tuple{Float64, Sampleable}}(
                    :high => (0.9, Delta(1.0)),
                    :low => (0.1, TruncatedNormal(0.05, 0.07, 0, 1))
                )
        else
            knockdown_dist = Dict{Symbol, Tuple{Float64, Sampleable}}(
                :high => (0.9, TruncatedNormal(0.90, 0.1, 0, 1)),
                :low => (0.1, TruncatedNormal(0.05, 0.07, 0, 1))
            )
        end

        lib = Library(knockdown_dist, max_phenotype_dists, kd_phenotype_relationships, crisprtype)

        for i in 1:num_runs
            push!(runs, (screen, lib, i))
        end
    end

    const overlap = intersect(fieldnames(FacsScreen), fieldnames(GrowthScreen))
    # custom function for handling both growth and a facs screen in the
    # same relational datastructure
    flatten_overlap = (setup, lib) -> begin
        local results = Any[typeof(setup)]
        for name in overlap
            push!(results, getfield(setup, name))
        end
        push!(results, typeof(setup) == GrowthScreen ? getfield(setup, :num_bottlenecks) : NaN)
        push!(results, as_array(lib)...)
        results
    end

    compare_reduction_methods = (bc_counts, genes) -> begin
        sig, noi = signal(bc_counts), noise(bc_counts)

        pvalmeanprod = auprc(abs.(genes[:pvalmeanprod_bin2_div_bin1]), genes[:class], Set([:increasing, :decreasing]))[1]
        pvalue = auprc(genes[:pvalue_bin2_div_bin1], genes[:class], Set([:increasing, :decreasing]))[1]
        effectsize = auprc(genes[:absmean_bin2_div_bin1], genes[:class], Set([:increasing, :decreasing]))[1]

        (sig, noi, pvalmeanprod, pvalue, effectsize)
    end

    before = time()
    results = pmap(args -> run_exp(args[1],
                                   args[2],
                                   compare_reduction_methods;
                                   flatten_func=flatten_overlap,
                                   run_idx=args[3]), runs)
    (!quiet) && println("$(time() - before) seconds")
    results = DataFrame(permutedims(hcat(results...), [2,1]))
    new_names = [:method; :score; :screentype; overlap...; :num_bottlenecks; array_names(Library)...; :run_idx]
    hierarchy = reshape([:signal, :noise, :product, :pvalue, :effectsize], 5, 1)
    results = construct_hierarchical_label(hierarchy, results, new_names)
    CSV.write(filepath, results)
end
