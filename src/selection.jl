"""
Given `cells`, a vector of integers, and `guides`, a vector of barcodes
performs simulated facs sorting of the cells into `bins` with the given
cutoffs. `σ` refers to the spread of the phenotype during the FACS screen.
"""
function facs_sort(cells::Vector{Int64}, guides::Vector{Barcode},
                   bins::Dict{Symbol, Tuple{Float64, Float64}}, σ::Float64)

    binned_cells = Array{Int64}[Int64[] for _ in length(keys(bins))]
    n_cells = length(cells)
    observed = zeros(n_cells)
    indices = collect(1:n_cells)
    for i in 1:n_cells
        cell = cells[i]
        observed[i] = rand(Normal(guides[cell].theo_phenotype, σ))
        guides[cell].obs_phenotype = observed[i]
    end
    sortperm!(indices, observed)
    permute!(cells, indices)
    results = Dict{Symbol, Vector{Int64}}()
    for (binname, cutoffs) in bins
        left = clamp(round(Int64, cutoffs[1]*n_cells), 1, n_cells)
        right = clamp(round(Int64, cutoffs[2]*n_cells), 1, n_cells)
        results[binname] = cells[left:right]
    end
    results
end
