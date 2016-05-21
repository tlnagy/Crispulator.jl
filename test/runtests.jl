using Base.Test
using DataStructures

include("../src/run.jl")

linear_response = response(Linear())

@test linear_response(0.0, 1.0) == 0.0
@test linear_response(1.0, 1.0) == 1.0
@test linear_response(0.0, -1.0) == 0.0
@test linear_response(1.0, -1.0) == -1.0

sigmoidal_response = response(Sigmoidal())

@test sigmoidal_response(0.0, 1.0) < 0.001
@test sigmoidal_response(1.0, 1.0) > 0.999
@test sigmoidal_response(0.0, -1.0) > -0.001
@test sigmoidal_response(1.0, -1.0) < -0.999

# Test the two different CRISPR types
lib = Library(CRISPRi())
setup = FacsScreen()
guides, guide_freqs_dist = setup_screen(setup, lib)
cells, cell_phenotypes = build_cells(lib.cas9_behavior, guides, guide_freqs_dist, 10^7)

function convert_cells_to_pop(cells, cell_phenotypes)
    cells_to_phenotypes = [DefaultDict(Float64, Int64, 0) for _ in 1:length(guides)]

    @inbounds for i in eachindex(cells)
        cells_to_phenotypes[cells[i]][cell_phenotypes[i]] += 1
    end

    cells_to_phenotypes
end

# In a CRISPRi screen all cells should have the same phenotype
@test all(map(length, convert_cells_to_pop(cells, cell_phenotypes)) .== 1)

lib = Library(CRISPRKO())
guides, guide_freqs_dist = setup_screen(setup, lib)
cells, cell_phenotypes = build_cells(lib.cas9_behavior, guides, guide_freqs_dist, 10^7)

arr = convert_cells_to_pop(cells, cell_phenotypes);
nonzeros = arr[find(x->length(x) != 1, arr)]
results = zeros(CRISPRKO().knockout_dist.K, length(nonzeros))
for (idx, guide) in enumerate(nonzeros)
    if -0.0 in keys(guide)
        @assert guide[0.0] == 0
        delete!(guide, 0.0)
        ks = sort(collect(keys(guide)))
        tot = sum(values(guide))
        results[:, idx] = [guide[key]/tot for key in ks]
    else
        @assert guide[0.0] != 0
        ks = sort(collect(keys(guide)), rev=true)
        tot = sum(values(guide))
        results[:, idx] = [guide[key]/tot for key in ks]
    end
end

@test abs(mean(results[1, :] ./ results[2, :]) - 1) < 0.025
@test abs(mean(results[1, :] ./ results[3, :])/4 - 1) < 0.025

println("All tests pass")
