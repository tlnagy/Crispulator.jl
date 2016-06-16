"""
Given `cells`, a vector of integers, and `guides`, a vector of barcodes
performs simulated facs sorting of the cells into `bins` with the given
cutoffs. `σ` refers to the spread of the phenotype during the FACS screen.
"""
function select(setup::FacsScreen,
                cells::Vector{Int64},
                cell_phenotypes::Vector{Float64},
                guides::Vector{Barcode};
                debug = false)

    σ = setup.σ
    bins = setup.bin_info
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

"""
Growth Screen selection
"""
function select(setup::GrowthScreen,
                initial_cells::AbstractArray{Int64},
                initial_cell_phenotypes::AbstractArray{Float64},
                guides::Vector{Barcode};
                debug=false)

    bottleneck_representation = setup.bottleneck_representation
    num_bottlenecks = setup.num_bottlenecks
    # all cells at all timepoints
    cellmat = zeros(Int64, length(guides)*bottleneck_representation)
    cpmat = zeros(Float64, size(cellmat))
    output_c = Array(Int64, length(initial_cells)*4);
    output_p = Array(Float64, size(output_c));
    cells = initial_cells # 1st timepoint slice
    picked = Array(Int64, size(cellmat, 1))
    cell_phenotypes = initial_cell_phenotypes

    for k in 1:num_bottlenecks
        num_inserted = grow!(cells, cell_phenotypes, output_c, output_p)
        cells, cell_phenotypes = sub(cellmat, :), sub(cpmat, :)
        sample!(1:num_inserted, picked, replace=false)
        copy!(sub(cellmat, :), sub(output_c, picked))
        copy!(sub(cpmat, :), sub(output_p, picked))
    end
    return Dict(:bin1 => initial_cells, :bin2 => cellmat)
end
