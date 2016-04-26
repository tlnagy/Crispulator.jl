"""
Constructs the guide library for `N` genes with `coverage` number of guides per
gene. Returns a tuple of guides and their relative frequencies (assigned randomly).
"""
function construct_library(N::Int64, coverage::Int64)
    # Assuming a high quality library has mostly good guides with some bad ones
    knockdown_model = MixtureModel([TruncatedNormal(0.97, 0.05, 0, 1), 
                                    TruncatedNormal(0.05, 0.07, 0, 1)], [0.9, 0.1])
    
    # The maximum expected phenotype is drawn from a mixture of three models.
    # we expect most genes to not have a phenotype, which is represented by a 
    # delta-like function centered a 0. 
    max_phenotype = MixtureModel([Normal(0, eps()), 
                                  TruncatedNormal(0.8, 0.2, 0.1, 1), 
                                  TruncatedNormal(-0.8, 0.2, -1, -0.1)], [0.8, 0.1, 0.1])
    # for the slopes we have equal chance for sigmoidal and linear options
    slope_model = MixtureModel([TruncatedNormal(2, 1.5, 0, 10), 
                                TruncatedNormal(30, 5, 10, Inf)])
    # inflections in the sigmoids occurs at high knockdowns
    inflection_model = TruncatedNormal(0.8, 0.2, 0, 1)
    
    barcodes = Barcode[]
    graphs = []
    
    for gene in 1:N
        
        # Generate knockdown-phenotype relationship for this gene
        params = [rand(max_phenotype),rand(slope_model),rand(inflection_model)]
        gene_response(x) = sigmoid(x, params...)
        
        for i in 1:coverage
            barcode_knockdown = rand(knockdown_model)
            # phenotype of barcode given its knockdown efficiency
            barcode_phenotype = gene_response(barcode_knockdown)
            (-0.01 < barcode_phenotype < 0.01) && (barcode_phenotype = 0.0)
            behavior = params[2] < 10 ? :linear : :sigmoidal
            if params[1] < -0.01
                class = :down
            elseif params[1] > 0.01
                class = :up
            else
                class = :inactive
            end
            push!(barcodes, Barcode(gene, barcode_knockdown, barcode_phenotype, behavior, class))
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
Given value(s) `x` apply the sigmoidal function with maximum
value `l`, a steepness `k`, and an inflection point `p`.
"""
sigmoid(x, l, k, p) = clamp(l./(1 + e.^(-k.*(x - p))), min(0, l), max(0, l))
