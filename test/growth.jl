using StatsBase

function test_grow_function(phenotypes, expected)
    @assert length(phenotypes) == length(expected)
    cells = collect(1:length(phenotypes))
    cell_phenotypes = phenotypes
    output_c = Array{Int}(length(cells)*4)
    output_p = Array{Float64}(size(output_c))

    num_runs = 10000
    data = Array{Int}(length(cells), num_runs)
    setup = GrowthScreen()
    setup.noise = 0.00001

    for i in 1:num_runs
        num_inserted = grow!(cells, cell_phenotypes, output_c, output_p, setup)
        data[:, i] = StatsBase.counts(output_c[1:num_inserted])
    end
    all(round.(mean(data, 2), 1) .== expected)
end

# test that, averaged over many runs, each phenotype corresponds to the
# expected number of doublings
@test test_grow_function([0, -1, 1, 0.5, -0.5], [2.0, 1.0, 4.0, 3.0, 1.5])
# unbalanced
@test test_grow_function([0, -1, -0.5, -0.8], [2.0, 1.0, 1.5, 1.2])
