
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
    for i in 1:length(cells)
        ρ::Float64 = cell_phenotypes[i]
        ρ_noisy = ρ + rand(noise_dist)
        decision = abs(ρ_noisy) < rand() ? 2 : 2^trunc(Int, 1 + sign(ρ_noisy))
        rng = num_inserted+1:num_inserted+decision
        output_c[rng] .= cells[i]
        output_p[rng] .= ρ
        num_inserted+=decision
    end
    num_inserted
end

"""
$(SIGNATURES)

Given a population `initial_cells` with `initial_cell_phenotypes`, this function
models the selection step of a [GrowthScreen](@ref). Cells are allowed to double
according to their growth phenotype and then the population is subsampled to
simulate passaging to maintain the `bottleneck_representation` as defined in
`setup`. The number of doublings + subsamplings are controlled by the
`num_bottlenecks` parameter in `setup`.

So for a screen where there are three bottlenecks and the
`length(initial_cells) == bottleneck_representation` then we would expect the
following behavior:

               Cell population
  ┌────────────────────────────────────────┐
2 │           .-           .r|          ..|│
  │         r*`|         _-' |        .r/ |│
  │       ."   |       .*    |      .r'   |│
  │    .-"      |    ,"`     |    ,/'     |│
  │  _-`        | .-/        |  .*`       |│
1 │_*           Lr'          |r/          |│
  └────────────────────────────────────────┘
   0                                      3
                 Bottleneck #

If the coverage of the library by `initial_cells` is less than the
`bottleneck_representation` than the population will be allowed to double until
passaging is needed or the `num_bottlenecks` is hit, which ever comes first.
"""
function select(setup::GrowthScreen,
                initial_cells::AbstractArray{Int},
                initial_cell_phenotypes::AbstractArray{Float64},
                guides::Vector{Barcode})

    bottleneck_representation = setup.bottleneck_representation
    num_bottlenecks = setup.num_bottlenecks

    # pre-allocated vectors
    cellmat = zeros(Int, length(guides)*bottleneck_representation)
    cpmat = zeros(Float64, size(cellmat))
    output_c = Array{Int}(undef, length(initial_cells)*4);
    output_p = Array{Float64}(undef, size(output_c));

    # 1st timepoint slice
    cells = initial_cells
    cell_phenotypes = initial_cell_phenotypes

    num_cells = 0
    picked = Array{Int}(undef, size(cellmat, 1))

    @debug "Initial cell population: $(length(initial_cells))"

    for k in 1:num_bottlenecks
        num_cells = grow!(cells, cell_phenotypes, output_c, output_p, setup)
        @debug "Doubling $k: cell population grew to $num_cells"

        # if we have more cells than we are planning to keep, select a random
        # subpopulation without replacement
        if num_cells > length(picked)
            sample!(1:num_cells, picked, replace=false)
            @debug "Selecting $(length(picked))"
        else
            # if we have room to grow, don't subsample, just shuffle and grow
            # the available space
            picked .= 0
            picked[1:num_cells] .= shuffle(1:num_cells)

            resize!(output_c, num_cells * 4)
            resize!(output_p, num_cells * 4)
        end

        ids = view(picked, picked .> 0)
        num_cells = length(ids)
        copy!(view(cellmat, 1:num_cells), view(output_c, ids))
        copy!(view(cpmat, 1:num_cells), view(output_p, ids))
        cells, cell_phenotypes = view(cellmat, 1:length(ids)), view(cpmat, 1:length(ids))
    end
    return OrderedDict(:bin1 => initial_cells, :bin2 => cellmat[1:num_cells])
end
