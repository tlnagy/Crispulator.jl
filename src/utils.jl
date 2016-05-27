function count(labels::AbstractArray{Symbol}, pos_labels::Set{Symbol})
    num_pos, num_neg = 0, 0
    for label in labels
        if label in pos_labels
            num_pos += 1
        else
            num_neg += 1
        end
    end
    num_pos, num_neg
end
"""
auroc(scores::AbstractArray{Float64}, classes::AbstractArray{Symbol}, pos_labels::Set{Symbol})

Optimized function for computing the area under the receiver operator characteristic
curve.
"""
function auroc(scores::AbstractArray{Float64}, classes::AbstractArray{Symbol}, pos_labels::Set{Symbol})
    num_scores = length(scores)
    ordering = sortperm(scores, rev=true)
    labels = classes[ordering]
    num_pos, num_neg = count(labels, pos_labels)

    tprs = Array(Float64, num_scores)
    fprs = Array(Float64, num_scores)

    auroc = 0.0
    tp = labels[1] in pos_labels ? 1 : 0
    fp = 1-tp
    tn, fn = num_neg - fp, num_pos - tp
    tprs[1] = tp/(tp+fn)
    fprs[1] = fp/(fp+tn)

    for i in 2:num_scores
        dtp = labels[i] in pos_labels ? 1 : 0
        tp += dtp
        fp += 1-dtp
        tn = num_neg - fp
        fn = num_pos - tp
        tprs[i] = tp/(tp+fn)
        fprs[i] = fp/(fp+tn)
        auroc += (tprs[i] + tprs[i-1])/2*(fprs[i] - fprs[i-1])
    end
    auroc, tprs, fprs
end

"""
auprc(scores::AbstractArray{Float64}, classes::AbstractArray{Symbol}, pos_labels::Set{Symbol})

Computes the area under the Precision-Recall curve using a lower
trapezoidal estimator, which is more accurate for skewed datasets.


K. Boyd, K. H. Eng, and C. D. Page, “Area under the Precision-Recall
Curve: Point Estimates and Confidence Intervals,” in Machine Learning
and Knowledge Discovery in Databases, H. Blockeel, K. Kersting,
S. Nijssen, and F. Železný, Eds. Springer Berlin Heidelberg, 2013,
pp. 451–466.
"""
function auprc(scores::AbstractArray{Float64}, classes::AbstractArray{Symbol}, pos_labels::Set{Symbol})
    num_scores = length(scores)
    ordering = sortperm(scores, rev=true)
    labels = classes[ordering]
    num_pos, num_neg = count(labels, pos_labels)

    tn, fn, tp, fp = 0, 0, num_pos, num_neg

    p = Array(Float64, num_scores)
    r = Array(Float64, num_scores)
    p[num_scores] = tp/(tp+fp)
    r[num_scores] = tp/(tp+fn)
    auprc, prev_r = 0.0, r[num_scores]
    pmin, pmax = p[num_scores], p[num_scores]

    # traverse scores from lowest to highest
    for i in num_scores-1:-1:1
        dtn = labels[i+1] in pos_labels ? 0 : 1
        tn += dtn
        fn += 1-dtn
        tp = num_pos - fn
        fp = num_neg - tn
        p[i] = tp/(tp+fp)
        r[i] = tp/(tp+fn)

        # update max precision observed for current recall value
        if r[i] == prev_r
            pmax = p[i]
        else
            pmin = p[i] # min precision is always at recall switch
            auprc += (pmin + pmax)/2*(prev_r - r[i])
            prev_r = r[i]
        end
    end
    auprc, p, r
end

"""
venn(scores::AbstractArray{Float64}, classes::AbstractArray{Symbol}, pos_labels::Set{Symbol})

Given `N` positive examples, computes the percentage of the top `N/2`
hits that are correct
"""
function venn(scores::AbstractArray{Float64}, classes::AbstractArray{Symbol}, pos_labels::Set{Symbol})
    num_scores = length(scores)
    ordering = sortperm(scores, rev=true)
    labels = classes[ordering]
    num_pos, num_neg = count(labels, pos_labels)

    n_venn = round(Int, num_pos/2)
    n_correct = 0

    for i in 1:n_venn
        (labels[i] in pos_labels) && (n_correct += 1)
    end
    n_correct/n_venn
end
