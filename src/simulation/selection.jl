import Base.Sort: select

"""
$(SIGNATURES)

Given `cells`, a vector of integers, and `guides`, a vector of barcodes,
performs simulated FACS sorting of the cells into `bins` with the given
cutoffs. `σ` refers to the spread of the phenotype during the FACS screen.
"""
function select(setup::FacsScreen,
                cells::Vector{Int},
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
    results = OrderedDict{Symbol, Vector{Int}}()

    for (binname, cutoffs) in bins
        left = clamp(round(Int, cutoffs[1]*n_cells), 1, n_cells)
        right = clamp(round(Int, cutoffs[2]*n_cells), 1, n_cells)
        results[binname] = cells[left:right]
    end
    results
end

function grow!(cells::AbstractArray{Int},
               cell_phenotypes::AbstractArray{Float64},
               output_c::AbstractArray{Int},
               output_p::AbstractArray{Float64},
               setup::GrowthScreen)
    num_inserted::Int = 0
    noise_dist = Normal(0, setup.noise)
    @inbounds for i in 1:length(cells)
        ρ::Float64 = cell_phenotypes[i]
        ρ_noisy = ρ + rand(noise_dist)
        decision = abs(ρ_noisy) < rand() ? 2 : 2^trunc(Int, 1 + sign(ρ_noisy))
        rng = num_inserted+1:num_inserted+decision
        output_c[rng] = cells[i]
        output_p[rng] = ρ
        num_inserted+=decision
    end
    num_inserted
end

"""
$(SIGNATURES)

Growth Screen selection
"""
function select(setup::GrowthScreen,
                initial_cells::AbstractArray{Int},
                initial_cell_phenotypes::AbstractArray{Float64},
                guides::Vector{Barcode};
                debug=false)

    bottleneck_representation = setup.bottleneck_representation
    num_bottlenecks = setup.num_bottlenecks
    # all cells at all timepoints
    cellmat = zeros(Int, length(guides)*bottleneck_representation)
    cpmat = zeros(Float64, size(cellmat))
    output_c = Array{Int}(length(initial_cells)*4);
    output_p = Array{Float64}(size(output_c));
    cells = initial_cells # 1st timepoint slice
    picked = Array{Int}(size(cellmat, 1))
    cell_phenotypes = initial_cell_phenotypes

    for k in 1:num_bottlenecks
        num_inserted = grow!(cells, cell_phenotypes, output_c, output_p, setup)
        cells, cell_phenotypes = view(cellmat, :), view(cpmat, :)
        sample!(1:num_inserted, picked, replace=false)
        copy!(view(cellmat, :), view(output_c, picked))
        copy!(view(cpmat, :), view(output_p, picked))
    end
    return OrderedDict(:bin1 => initial_cells, :bin2 => cellmat)
end
