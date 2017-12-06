# test AUROC compute functions
tprs = zeros(101)
fprs = zeros(101)

scores = collect(1:-0.01:0)
classes = [:b,:a,:a,:a,:b,:a,:a,:b,:b,:b,:a,:a,:b,:b,:b,:b,:a,:b,:a,:a,
    :a,:b,:a,:a,:b,:a,:b,:a,:a,:a,:b,:b,:a,:a,:a,:a,:a,:b,:b,:b,:a,:b,
    :a,:a,:b,:a,:a,:a,:a,:b,:a,:a,:b,:a,:a,:a,:a,:a,:b,:a,:b,:b,:a,:b,
    :a,:a,:a,:b,:b,:a,:a,:b,:b,:a,:b,:b,:b,:a,:b,:a,:b,:b,:a,:a,:b,:b,
    :b,:a,:b,:b,:b,:a,:a,:a,:a,:a,:a,:b,:a,:a,:b]
@test isapprox(auroc(scores, classes, Set([:a]))[1], 0.513556, atol=1e-6)
@test venn(scores, classes, Set([:a])) == 0.5357142857142857

scores = collect(1:-0.1:0)
classes = [:b, :b, :b, :b, :a, :a, :a, :a, :a, :a, :a]
@test auroc(scores, classes, Set([:a]))[1] == 0.0
@test venn(scores, classes, Set([:a])) == 0.0

scores = collect(1:-0.1:0)
classes = [:a, :a, :a, :a, :a, :a, :a, :a, :b, :b, :b]
@test auroc(scores, classes, Set([:a]))[1] == 1.0
@test venn(scores, classes, Set([:a])) == 1.0

# this algorithm is a direct implementation from the Boyd paper
function compute_auprc_posthoc(p, r)
    unqs = unique(r)
    pmins = p[[findlast(r, unq) for unq in unqs]]
    pmaxes = p[[findfirst(r, unq) for unq in unqs]]
    n = length(unqs)
    auprc = 0.0
    for i in n-1:-1:1
        auprc += (pmins[i] + pmaxes[i+1])/2 * (unqs[i+1] - unqs[i])
    end
    auprc
end

# Test multiple positive labels
classes = [:c,:a,:a,:b,:b,:b,:c,:c,:a,:c,:a,:a]
scores = Float64[80,5,157,169,158,166,115,18,6,78,31,5]
@test venn(scores, classes, Set([:a])) == 0.0
@test venn(scores, classes, Set([:b])) == 1.0
@test venn(scores, classes, Set([:a, :b])) == 1.0

@test auroc(scores, classes, Set([:a]))[1] == 0.14285714285714288
@test auroc(scores, classes, Set([:b]))[1] == 1.0
@test auroc(scores, classes, Set([:a, :b]))[1] == 0.53125

@test auprc(scores, classes, Set([:b]))[1] == 1.0
a, p, r = auprc(scores, classes, Set([:a, :b]))
@test a == compute_auprc_posthoc(p, r)
a, p, r = auprc(scores, classes, Set([:a]))
@test a == compute_auprc_posthoc(p, r)

# test for #30
classes = [:a, :a, :a, :b, :c, :a, :c, :c, :c, :b, :b]
scores = [-10., -8, -6, -4, -2, -1, 1, 1, 2, 3, 7]

@test venn(scores, classes, Set([:a]), rev=true) == 0.0
@test venn(scores, classes, Set([:a]), rev=false) == 1.0
@test venn(abs.(scores), classes, Set([:a, :b]), rev=true) == 1.0

@test auprc(scores, classes, Set([:a]), rev=true)[1] ≈ 0.21246843434343435
@test auprc(scores, classes, Set([:a]), rev=false)[1] ≈ 0.9083333333333333
@test auprc(abs.(scores), classes, Set([:a, :b]), rev=true)[1] ≈ 0.9662698412698412

function compute_bias(recalls, precisions, X, Y)
    true_auprc = 0.0
    for i in 2:length(recalls)
        true_auprc += (precisions[i-1] + precisions[i])/2*(recalls[i - 1] - recalls[i])
    end

    num_runs = 1000
    auprc_sum = 0.0
    cat_dist = Categorical([1-π, π])
    @inbounds for run in 1:num_runs
        cat = rand(cat_dist, 10^3)
        x, y = StatsBase.counts(cat, 1:2)
        scores = [rand(X, x); rand(Y, y)]
        auprc_sum += auprc(scores, [repmat(classes[1:1], x); repmat(classes[2:2], y)], Set([:b]))[1]
    end
    isapprox(auprc_sum/num_runs, true_auprc, atol=0.025)
end

π = 0.1
test_dists = Array[
     [Normal(0, 1), Normal(1,1)],
     [Beta(2, 5), Beta(5, 2)],
     [Uniform(0, 1), Uniform(0.5, 1.5)]
]
# true precision, recall functions
recall(xs, Y) = 1-cdf.(Y, xs)
function precision(xs, π, X, Y)
    res = π*recall(xs, Y)./(π*recall(xs, Y) + (1-π)*(1-cdf.(X, xs)))
    res[isnan.(res)] = 1
    res
end
xs = linspace(-10, 10, 1000)

for dists in test_dists
    X, Y = dists
    recalls = recall(xs, Y)
    precisions = precision(xs, π, X, Y)
    classes = [:a, :b]
    @test compute_bias(recalls, precisions, X, Y)
end
