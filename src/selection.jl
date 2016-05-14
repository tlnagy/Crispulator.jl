"""
Given `cells`, a vector of integers, and `guides`, a vector of barcodes
performs simulated facs sorting of the cells into `bins` with the given
cutoffs. `σ` refers to the spread of the phenotype during the FACS screen.
"""
function facs_sort(cells::Vector{Int64}, guides::Vector{Barcode},
                   bins::Dict{Symbol, Tuple{Float64, Float64}}, σ::Float64)

    n_cells = length(cells)
    observed = zeros(n_cells)
    for i in 1:n_cells
        cell = cells[i]
        observed[i] = rand(Normal(guides[cell].theo_phenotype, σ))
        guides[cell].obs_phenotype = observed[i]
    end
    indices = sortperm(observed)
    cells = cells[indices]
    results = Dict{Symbol, Vector{Int64}}()

    for (binname, cutoffs) in bins
        left = clamp(round(Int64, cutoffs[1]*n_cells), 1, n_cells)
        right = clamp(round(Int64, cutoffs[2]*n_cells), 1, n_cells)
        results[binname] = cells[left:right]
    end
    results
end

function grow!(cells::AbstractArray{Int64}, guides::Vector{Barcode}, output)
    num_inserted::Int = 0
    @inbounds for i in 1:length(cells)
        id::Int64 = cells[i]
        ρ::Float64 = guides[id].theo_phenotype
        decision = abs(ρ) < rand() ? 2 : 2^trunc(Int, 1 + sign(ρ))
        output[num_inserted+1:num_inserted+decision] = id
        num_inserted+=decision
    end
    num_inserted
end

function growth_assay(initial_cells::AbstractArray{Int64},
                      guides::Vector{Barcode},
                      num_bottlenecks::Int64,
                      bottleneck_representation::Int64;
                      debug=false)

    # all cells at all timepoints
    cellmat = zeros(Int64, length(guides)*bottleneck_representation, num_bottlenecks)
    output = Array(Int64, length(initial_cells)*4);
    cells = initial_cells # 1st timepoint slice

    for k in 1:num_bottlenecks
        num_inserted = grow!(cells, guides, output)
        cells = sub(cellmat, :, k)
        sample!(sub(output, 1:num_inserted), cells, replace=false)
    end
    if debug
        d = Dict([symbol("bin", i+1)=>cellmat[:, i] for i in 1:num_bottlenecks])
        d[:bin1] = initial_cells
        return d
    else
        return Dict(:bin1 => initial_cells, :bin2 => cellmat[:, num_bottlenecks])
    end
end
