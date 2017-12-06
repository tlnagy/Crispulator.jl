function testselection(setup)
    N = setup.num_genes
    guides = Barcode[]

    for (gene, categ) in zip(1:N, [:decreasing, :inactive, :increasing])
        push!(guides, Barcode(gene, 1, gene-2, :linear, categ))
    end

    expand_to = setup.bottleneck_representation * length(guides)

    initial_cells = Array{Int}(expand_to)
    cell_phenotypes = Array{Float64}(size(initial_cells))
    for i in 1:length(initial_cells)
        initial_cells[i] = mod1(i, N)
        cell_phenotypes[i] = guides[initial_cells[i]].theo_phenotype
    end

    # ensure guide level phenotype isn't used so that CRISPRn works
    for guide in guides
        guide.theo_phenotype = -Inf
    end

    num_runs = 100
    results = Array{Float64}(N, num_runs)

    for i in 1:num_runs
        data = select(setup, initial_cells, cell_phenotypes, guides)
        data1 = StatsBase.counts(data[:bin1], 1:N)
        data2 = StatsBase.counts(data[:bin2], 1:N)
        results[:, i] = log2.(data2 ./ data1)
    end
    results
end

function testfacs(σ, quant, width)
    N = 3
    setup = FacsScreen()
    setup.num_genes = N
    setup.coverage = 1
    setup.bottleneck_representation = 20000
    setup.bin_info = Dict{Symbol, Tuple{Float64, Float64}}(:bin1 => (0, width), :bin2 => (1-width, 1))
    setup.σ = σ

    results = testselection(setup)

    ideal = log2(cdf.(Normal(-1, σ), quant)/cdf.(Normal(1,σ), quant))
    tocompare = collect(zip(mean(results, 2), [-ideal, 0.0, ideal]))
    all(Bool[isapprox(x[1], x[2], atol=0.1) for x in tocompare])
end

function testgrowth(ideal, num)
    N = 3
    setup = GrowthScreen()
    setup.num_genes = N
    setup.coverage = 1
    setup.num_bottlenecks = num
    setup.bottleneck_representation = 1000

    results = testselection(setup)
    tocompare = collect(zip(mean(results, 2), log2.(ideal)))
    all(Bool[isapprox(x[1], x[2], atol=0.1) for x in tocompare])
end

begin
    srand(1)
    @test testfacs(0.5, -0.50138209, 1/3)
    @test testfacs(0.5, 0, 1/2)
    @test testfacs(1, -0.58093999, 1/3)

    @test testgrowth([1/7, 2/7, 4/7]./(1/3), 1)
    @test testgrowth([1/21, 4/21,16/21]./(1/3), 2)
    srand()
end
