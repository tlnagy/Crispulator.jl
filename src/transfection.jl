function transfect(guides::Vector{Barcode}, guide_freqs_dist::Distributions.Categorical,
                   cell_count::Int64, moi::Float64, expand_to::Int64)
    cells = rand(guide_freqs_dist, round(Int64, pdf(Poisson(moi), 1)*cell_count))
    num_cells = length(cells)

    multiples = 1
    if expand_to > num_cells
        multiples = ceil(Int64, expand_to/num_cells)
        expansion = Array(Int64, num_cells*multiples)
        for rep in 1:multiples
            @inbounds expansion[(rep-1)*num_cells+1:rep*num_cells] = cells
        end
        cells = expansion
    else
        warn("targeted expansion is smaller than initial cell count" *
             ", using $num_cells instead of $expand_to")
    end

    initial_freqs = StatsBase.counts(cells, 1:length(guides)) ./ length(cells)

    for i in 1:length(guides)
        @inbounds guides[i].initial_freq = initial_freqs[i]
    end
    cells
end
