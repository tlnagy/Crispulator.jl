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
        cells = sample(cells, expand_to)
    end

    initial_freqs = StatsBase.counts(cells, 1:length(guides)) ./ length(cells)

    for i in 1:length(guides)
        @inbounds guides[i].initial_freq = initial_freqs[i]
    end
    cells
end

function transfect(setup::GrowthScreen,
                   guides::Vector{Barcode},
                   guide_freqs_dist::Categorical)

    num_guides = length(guides)
    cell_count = num_guides * setup.representation
    initial_cells = rand(guide_freqs_dist, round(Int64, pdf(Poisson(setup.moi), 1)*cell_count))
    target = num_guides * setup.bottleneck_representation

    if target < length(initial_cells)
        cells = sample(initial_cells, target)
        num_doublings = -1
    else
        cells = copy(initial_cells)
        num_inserted = length(cells)
        num_doublings = 0
        output = Array(Int64, target*4);

        while num_inserted < target
            num_inserted = grow!(cells, guides, output)
            cells = copy(sub(output, 1:num_inserted))
            num_doublings += 1
        end
    end

    initial_freqs = counts(cells, 1:length(guides)) ./ length(cells)

    for i in 1:length(guides)
        @inbounds guides[i].initial_freq = initial_freqs[i]
    end
    cells, num_doublings
end
