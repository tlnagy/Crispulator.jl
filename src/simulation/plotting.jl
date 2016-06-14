using Gadfly

function lib_plot(lib::Library)
    draw(PNG("plots/library/guide_knockdown.png", 12cm, 10cm, dpi=300),
    plot(x=rand(lib.knockdown_dist, 10^6), Geom.histogram,
    Guide.xlabel("Gene knockdown"),
    Guide.title("Library-wide guide knockdown efficiencies")))

    num_steps = 100
    xs_range = linspace(0, 1, num_steps)
    max_rep = 10
    xs = Array(Float64, num_steps*max_rep)
    ys = Array(Float64, length(xs))
    reps = Array(Int64, length(xs))
    behaviors = Array(Symbol, length(xs))

    classes, max_phenotypes = [], []
    rep = 1
    for gene in [rand_gene(lib) for _ in 1:10^4]
        push!(classes, gene[1])
        push!(max_phenotypes, gene[2])
        if rep <= max_rep
            xs[(rep-1)*num_steps+1:rep*num_steps] = xs_range
            gene_response = response(gene[4])
            ys[(rep-1)*num_steps+1:rep*num_steps] = gene_response(collect(xs_range), gene[2])
            behaviors[(rep-1)*num_steps+1:rep*num_steps] = gene[3]
            reps[(rep-1)*num_steps+1:rep*num_steps] = rep
        end
        rep += 1
    end

    result = DataFrame(class=classes, max_phenotype=max_phenotypes)
    draw(PNG("plots/library/max_phenotypes.png", 12cm, 10cm, dpi=300),
    plot(result, x=:max_phenotype, color=:class, Geom.histogram,
    Guide.xlabel("Theoretical Phenotype"),
    Guide.title("Distribution of phenotypes within response categories")))

    result = DataFrame(xs=xs, ys=ys, behaviors=behaviors, reps=reps)

    draw(PNG("plots/library/kd_phenotype_relationships.png", 12cm, 10cm, dpi=300),
    plot(result, x=:xs, y=:ys, color=:behaviors, Geom.line))
end
