# test AUROC compute functions
tprs = zeros(101)
fprs = zeros(101)

scores = collect(1:-0.01:0)
classes = [:inactive,:a,:a,:a,:inactive,:a,:a,:inactive,:inactive,:inactive,
    :a,:a,:inactive,:inactive,:inactive,:inactive,:a,:inactive,:a,:a,:a,
    :inactive,:a,:a,:inactive,:a,:inactive,:a,:a,:a,:inactive,:inactive,:a,:a,
    :a,:a,:a,:inactive,:inactive,:inactive,:a,:inactive,:a,:a,:inactive,:a,:a,
    :a,:a,:inactive,:a,:a,:inactive,:a,:a,:a,:a,:a,:inactive,:a,:inactive,
    :inactive,:a,:inactive,:a,:a,:a,:inactive,:inactive,:a,:a,:inactive,
    :inactive,:a,:inactive,:inactive,:inactive,:a,:inactive,:a,:inactive,
    :inactive,:a,:a,:inactive,:inactive,:inactive,:a,:inactive,:inactive,
    :inactive,:a,:a,:a,:a,:a,:a,:inactive,:a,:a,:inactive]
@test isapprox(compute_auroc(scores, classes, 101), 0.513556, atol=1e-6)
@test isapprox(compute_roc!(scores, classes, 101, tprs, fprs), 0.513556, atol=1e-6)

scores = collect(1:-0.1:0)
classes = [:inactive, :inactive, :inactive, :inactive, :a, :a, :a, :a, :a, :a, :a]
@test compute_auroc(scores, classes, 10) == 0.0
@test compute_roc!(scores, classes, 10, tprs, fprs) == 0.0

scores = collect(1:-0.1:0)
classes = [:a, :a, :a, :a, :a, :a, :a, :a, :inactive, :inactive, :inactive]
@test compute_auroc(scores, classes, 10) == 1.0
@test compute_roc!(scores, classes, 10, tprs, fprs) == 1.0

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
        auprc_sum += fast_auprc(scores, [repmat(classes[1:1], x); repmat(classes[2:2], y)], :b)[1]
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
recall(xs, Y) = 1-cdf(Y, xs)
function precision(xs, π, X, Y)
    res = π*recall(xs, Y)./(π*recall(xs, Y) + (1-π)*(1-cdf(X, xs)))
    res[isnan(res)] = 1
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
