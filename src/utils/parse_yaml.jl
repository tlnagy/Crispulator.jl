function Base.parse(data::Dict{Any, Any})

    # library parameters
    libdata = data["library"]
    genome = libdata["genome"]
    guides = libdata["guides"]

    num_genes, coverage = -1, -1
    frac_inc_genes, frac_dec_genes = 0.0, 0.0
    crisprtype = "CRISPRn"
    frac_hq, mean_kd = 0.0, 0.0

    try
        num_genes = genome["num-genes"]::Int64
        coverage = genome["num-guides-per-gene"]::Int64
        frac_inc_genes = genome["frac-increasing-genes"]::Float64
        frac_dec_genes = genome["frac-decreasing-genes"]::Float64

        crisprtype = guides["crispr-type"]
        frac_hq = guides["frac-high-quality"]::Float64
        mean_kd = guides["mean-high-quality-kd"]::Float64

    catch e
        if isa(e, InexactError)
            error("Number of genes and coverage should both be integers")
        else
            rethrow(e)
        end
    end

    if lowercase(crisprtype) == "crispri"
        cas9_behavior = CRISPRi()
        knockdown_dist = Dict{Symbol, Tuple{Float64, Sampleable}}(
            :high => (frac_hq, TruncatedNormal(mean_kd, 0.1, 0, 1)),
            :low => (1-frac_hq, TruncatedNormal(0.05, 0.07, 0, 1))
        )
    elseif lowercase(crisprtype) == "crisprn"
        cas9_behavior = CRISPRKO()
        knockdown_dist = Dict{Symbol, Tuple{Float64, Sampleable}}(
            :high => (frac_hq, Delta(1.0)),
            :low => (1-frac_hq, TruncatedNormal(0.05, 0.07, 0, 1))
        )
    else
        error("CRISPR type must be either CRISPRi or CRISPRn")
    end

    max_phenotype_dists = Dict{Symbol, Tuple{Float64, Sampleable}}(
        :inactive => (1-frac_inc_genes-frac_dec_genes-0.05, Delta(0.0)),
        :negcontrol => (0.05, Delta(0.0)),
        :increasing => (frac_inc_genes, TruncatedNormal(0.1, 0.1, 0.025, 1)),
        :decreasing => (frac_dec_genes, TruncatedNormal(-0.55, 0.2, -1, -0.1))
    )
    lib = Library(max_phenotype_dists, knockdown_dist, cas9_behavior)

    println(data["screen"])

end
