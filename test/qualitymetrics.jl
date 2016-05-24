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
