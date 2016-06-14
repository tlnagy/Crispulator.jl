# There are many 3 primary stages of the experiment where representation
# is important: at transfection, at the bottleneck(s), and at
# sequencing.

function main(filepath; debug=false)

    if !debug
        parameters = Dict{Symbol, Vector}(
            :representation => map(x->round(Int64, x), logspace(0, 4, 20)),
            :bottleneck_representation => map(x->round(Int64, x),  logspace(0,4,20)),
            :seq_depth => [10^2, 10^3, 10^4]
        )
        num_runs = 25
    else
        parameters = Dict{Symbol, Vector}(
            :representation => map(x->round(Int64, x), logspace(0, 4, 2)),
            :bottleneck_representation => map(x->round(Int64, x),  logspace(0,4,2)),
            :seq_depth => [10^4]
        )
        num_runs = 1
    end

    const overlap = intersect(fieldnames(FacsScreen), fieldnames(GrowthScreen))
    # custom function for handling both growth and a facs screen in the
    # same relational datastructure
    flatten_overlap = (screen) -> begin
        local results = Any[typeof(screen)]
        for name in overlap
            push!(results, getfield(screen, name))
        end
        results
    end

    # computes the auprcs of increasing and decreasing genes separately
    get_auprcs = genes -> begin
        a = auprc(abs(genes[:pvalmeanprod]), genes[:class], Set([:increasing, :decreasing]))
        d = auprc(genes[:pvalmeanprod], genes[:class], Set([:decreasing]), rev=false)
        i = auprc(genes[:pvalmeanprod], genes[:class], Set([:increasing]))
        (a[1], d[1], i[1])
    end

    results = []
    for screentype in [FacsScreen(), GrowthScreen()]
        for crisprtype in [CRISPRi(), CRISPRKO()]
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

            result = @time pmap(args -> run_exp(args[1], lib, get_auprcs; run_idx=args[2], flatten_func=flatten_overlap), runs)
            result = DataFrame(hcat(result...)')
            result[:crisprtype] = typeof(crisprtype)
            push!(results, result)
        end
    end
    results = vcat(results...)
    names!(results, [:all; :inc; :dec; :screen; overlap...; :run; :crisprtype])

    writetable(filepath, results)
end
