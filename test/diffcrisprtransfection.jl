function convert_cells_to_pop(cells, cell_phenotypes, guides)
    cells_to_phenotypes = [DefaultDict{Float64, Int}(0) for _ in 1:length(guides)]

    @inbounds for i in eachindex(cells)
        cells_to_phenotypes[cells[i]][cell_phenotypes[i]] += 1
    end

    cells_to_phenotypes
end

function test_crispri_construction()
    # Test the two different CRISPR types
    lib = Library(CRISPRi())
    setup = FacsScreen()
    guides, guide_freqs_dist = construct_library(setup, lib)
    cells, cell_phenotypes = build_cells(lib.cas9_behavior, guides, guide_freqs_dist, 10^6)

    # In a CRISPRi screen all cells should have the same phenotype
    all(map(length, convert_cells_to_pop(cells, cell_phenotypes, guides)) .== 1)
end
@test test_crispri_construction()

function test_crisprko_construction()
    lib = Library(CRISPRn())
    setup = FacsScreen()
    guides, guide_freqs_dist = construct_library(setup, lib)
    cells, cell_phenotypes = build_cells(lib.cas9_behavior, guides, guide_freqs_dist, 10^7)

    arr = convert_cells_to_pop(cells, cell_phenotypes, guides)
    nonzeros = arr[find(x->length(x) != 1, arr)]
    results = zeros(CRISPRn().knockout_dist.K, length(nonzeros))
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

    isapprox(mean(results[1, :] ./ results[2, :]), 1, atol=0.25) &&
    isapprox(mean(results[1, :] ./ results[3, :]), 4, atol=0.25)
end

@test test_crisprko_construction()
