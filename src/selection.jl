"""
Given `cells`, a vector of integers, and `guides`, a vector of barcodes
performs simulated facs sorting of the cells into `bins` with the given
cutoffs. `σ` refers to the spread of the phenotype during the FACS screen.
"""
function facs_sort(cells::Vector{Int64}, cell_phenotypes::Vector{Float64},
                   guides::Vector{Barcode},
                   bins::Dict{Symbol, Tuple{Float64, Float64}}, σ::Float64)

    @assert size(cells) == size(cell_phenotypes)
    n_cells = length(cells)
    observed = zeros(n_cells)
    @inbounds for i in 1:n_cells
        observed[i] = rand(Normal(cell_phenotypes[i], σ))
        guides[cells[i]].obs_phenotype = observed[i]
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

function grow!(cells::AbstractArray{Int64}, cell_phenotypes::AbstractArray{Float64},
               output_c::AbstractArray{Int64}, output_p::AbstractArray{Float64})
    num_inserted::Int = 0
    @inbounds for i in 1:length(cells)
        ρ::Float64 = cell_phenotypes[i]
        decision = abs(ρ) < rand() ? 2 : 2^trunc(Int, 1 + sign(ρ))
        rng = num_inserted+1:num_inserted+decision
        output_c[rng] = cells[i]
        output_p[rng] = ρ
        num_inserted+=decision
    end
    num_inserted
end

function growth_assay(initial_cells::AbstractArray{Int64},
                      initial_cell_phenotypes::AbstractArray{Float64},
                      guides::Vector{Barcode},
                      num_bottlenecks::Int64,
                      bottleneck_representation::Int64;
                      debug=false)

    # all cells at all timepoints
    cellmat = zeros(Int64, length(guides)*bottleneck_representation, num_bottlenecks)
    cpmat = zeros(Float64, size(cellmat))
    output_c = Array(Int64, length(initial_cells)*4);
    output_p = Array(Float64, size(output_c));
    cells = initial_cells # 1st timepoint slice
    picked = Array(Int64, size(cellmat, 1))
    cell_phenotypes = initial_cell_phenotypes

    for k in 1:num_bottlenecks
        num_inserted = grow!(cells, cell_phenotypes, output_c, output_p)
        cells, cell_phenotypes = sub(cellmat, :, k), sub(cpmat, :, k)
        sample!(collect(1:num_inserted), picked, replace=false)
        copy!(sub(cellmat, :, k), output_c[picked])
        copy!(sub(cpmat, :, k), output_p[picked])
    end
    if debug
        d = Dict([symbol("bin", i+1)=>cellmat[:, i] for i in 1:num_bottlenecks])
        d[:bin1] = initial_cells
        return d
    else
        return Dict(:bin1 => initial_cells, :bin2 => cellmat[:, num_bottlenecks])
    end
end
