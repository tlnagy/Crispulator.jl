"""
Given value(s) `x` apply the sigmoidal function with maximum
value `l`, a steepness `k`, and an inflection point `p`.
"""
sigmoid(x, l, k, p) = clamp(l./(1 + e.^(-k.*(x - p))), min(0, l), max(0, l))

"""
Given value(s) `x` apply a simple linear function with maximum value `l`
"""
linear(x, l) = clamp(l.*x, min(0, l), max(0, l))

abstract KDPhenotypeRelationship

type Linear <: KDPhenotypeRelationship end

type Sigmoidal <: KDPhenotypeRelationship
    "Slopes of the sigmoid are pulled from this distribution"
    slope_dist::Distribution

    "Inflection points are pulled from this distribution"
    inflection_dist::Distribution

    function Sigmoidal()
        # slopes are pretty steep so that there is clear differentiation
        # between linear and sigmoidal
        slopes = TruncatedNormal(30, 5, 10, Inf)
        # inflections are centered at high KD
        inflections = TruncatedNormal(0.8, 0.2, 0, 1)
        new(slopes, inflections)
    end
end

response(::Linear) = linear
function response(sig::Sigmoidal)
    slope, inflection = rand(sig.slope_dist), rand(sig.inflection_dist)
    (x, l) -> sigmoid(x, l, slope, inflection)
end

"""
Wrapper containing all library construction parameters
"""
type Library
    "Distribution of guide knockdown efficiencies"
    knockdown_dist::Distribution

    """
    Maximum phenotype categories mapped to their probability of being
    selected and the distribution to draw from if they are selected
    """
    max_phenotype_dists::Dict{Int64, Tuple{Symbol, Sampleable}}
    phenotype_probs::Categorical

    """
    Knockdown-phenotype relationships mapped to their probability of
    being selected and their respective KDPhenotypeRelationship
    """
    kd_phenotype_relationships::Dict{Int64, Tuple{Symbol, KDPhenotypeRelationship}}
    relationship_probs::Categorical

    function Library(knockdown_dist::Distribution,
                     max_phenotype_dists::Dict{Symbol, Tuple{Float64, Sampleable}},
                     kd_phenotype_relationships::Dict{Symbol, Tuple{Float64, KDPhenotypeRelationship}})

        max_p = unroll(max_phenotype_dists)
        rela = unroll(kd_phenotype_relationships)
        new(knockdown_dist, max_p[1], max_p[2], rela[1], rela[2])
    end
end

function Library()
    max_phenotype_dists = Dict{Symbol, Tuple{Float64, Sampleable}}(
        :inactive => (0.75, Delta(0.0)),
        :negcontrol => (0.05, Delta(0.0)),
        :increasing => (0.1, TruncatedNormal(0.55, 0.2, 0.1, 1)),
        :decreasing => (0.1, TruncatedNormal(-0.55, 0.2, -1, -0.1))
    )
    Library(max_phenotype_dists)
end

function Library(max_phenotype_dists::Dict{Symbol, Tuple{Float64, Sampleable}})
    # Assuming a high quality library has mostly good guides with some bad ones
    knockdown_dist = MixtureModel([TruncatedNormal(0.90, 0.1, 0, 1),
                                   TruncatedNormal(0.05, 0.07, 0, 1)], [0.9, 0.1])

    kd_phenotype_relationships = Dict{Symbol, Tuple{Float64, KDPhenotypeRelationship}}(
        :linear => (0.75, Linear()),
        :sigmoidal => (0.25, Sigmoidal())
    )
    Library(knockdown_dist, max_phenotype_dists, kd_phenotype_relationships)
end

function unroll{T}(data::Dict{Symbol, Tuple{Float64, T}})
    probs = Float64[]
    results = Dict{Int64, Tuple{Symbol, T}}()
    count = 0
    for (key, value) in data
        push!(probs, value[1])
        results[count+=1] = (key, value[2])
    end
    results, Categorical(probs)
end

function rand_gene(lib::Library)
    class, dist = lib.max_phenotype_dists[rand(lib.phenotype_probs)]
    max_phenotype = rand(dist)
    behavior, relationship = lib.kd_phenotype_relationships[rand(lib.relationship_probs)]
    class, max_phenotype, behavior, relationship
end

"""
Constructs the guide library for `N` genes with `coverage` number of guides per
gene. Returns a tuple of guides and their relative frequencies (assigned randomly).
"""
function construct_library(lib::Library, N::Int64, coverage::Int64)
    barcodes = Barcode[]

    for gene in 1:N

        gene_class, max_phenotype, gene_behavior, relationship = rand_gene(lib)
        # get this gene's KD-phenotype relationship
        kd_response = response(relationship)

        for i in 1:coverage
            barcode_knockdown = rand(lib.knockdown_dist)
            # phenotype of barcode given its knockdown efficiency
            barcode_phenotype = kd_response(barcode_knockdown, max_phenotype)
            push!(barcodes, Barcode(gene, barcode_knockdown, barcode_phenotype,
                  gene_behavior, gene_class))
        end
    end

    # the expected z-score for a 90% confidence interval is 2x1.645=3.29
    # construct a normal distribution with a σ = 1/3.29 so that there is
    # a 10-fold difference between the 95th/5th percentiles
    vals = 10.^rand(Normal(0, 1/3.29), N*coverage)
    guide_freqs = vals/sum(vals)
    barcodes, guide_freqs
end

"""
Constructs a delta function at a given δ value
"""
type Delta <: Distributions.Sampleable{Univariate, Discrete}
    δ::Float64
end

Base.rand(d::Delta) = d.δ
