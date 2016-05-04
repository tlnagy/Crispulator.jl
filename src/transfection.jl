function transfect(guides::Vector{Barcode}, guide_freqs_dist::Distributions.Categorical,
                   cell_count::Int64, moi::Float64)
    cells = rand(guide_freqs_dist, round(Int64, pdf(Poisson(moi), 1)*cell_count))
    num_cells = length(cells)

    # simple linear expansion of 10-fold
    # Maybe play around with alternative strategies here
    expansion = Array(Int64, num_cells*10)
    for rep in 1:10
        @inbounds expansion[(rep-1)*num_cells+1:rep*num_cells] = cells
    end
    cells = expansion

    initial_freqs = StatsBase.counts(cells, 1:length(guides)) ./ length(cells)

    for i in 1:length(guides)
        @inbounds guides[i].initial_freq = initial_freqs[i]
    end
    cells
end
