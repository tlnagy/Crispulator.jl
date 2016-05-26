function compute_roc(df::DataFrame, score_column::Symbol, numsteps::Int64)
    tprs = zeros(numsteps)
    fprs = zeros(numsteps)
    auroc = compute_roc!(Vector(df[score_column]), Vector(df[:class]) , numsteps, tprs, fprs)
    return auroc, tprs, fprs
end

function _binary_classifier(scores, classes, threshold)
    tp, fp, tn, fn = 0,0,0,0

    @inbounds @fastmath for i in 1:length(scores)
        if scores[i] >= threshold
            if classes[i] != :inactive
                tp += 1
            else
                fp += 1
            end
        else
            if classes[i] == :inactive
                tn += 1
            else
                fn += 1
            end
        end
    end
    tp, fp, tn, fn
end

"""
Optimized function for computing the receiver operator characteristic curve.
"""
function compute_roc!(scores::Vector{Float64}, classes::Vector{Symbol}, numsteps::Int64,
                      tprs::Vector{Float64}, fprs::Vector{Float64})

    num_scores = length(scores)
    min_value, max_value = extrema(scores)
    thresholds = linspace(max_value, min_value, numsteps)
    auroc = 0.0

    @inbounds for (idx, threshold) in enumerate(thresholds)

        tp, fp, tn, fn = _binary_classifier(scores, classes, threshold)

        tprs[idx] = tp/(tp+fn)
        fprs[idx] = fp/(fp+tn)
        if idx > 1
            auroc += (tprs[idx] + tprs[idx-1])/2 * (fprs[idx] - fprs[idx-1])
        end
    end
    auroc
end

"""
Optimized function for computing the area under the receiver operator characteristic
curve.
"""
function compute_auroc(scores::Vector{Float64}, classes::Vector{Symbol}, numsteps::Int64)

    num_scores, idx = length(scores), 1
    min_value, max_value = extrema(scores)
    thresholds = linspace(max_value, min_value, numsteps)
    auroc, prev_tpr, prev_fpr = 0.0, 0.0, 0.0

    @inbounds for threshold in thresholds

        tp, fp, tn, fn = _binary_classifier(scores, classes, threshold)

        tpr = tp/(tp+fn)
        fpr = fp/(fp+tn)
        if idx > 1
            auroc += (tpr + prev_tpr)/2 * (fpr - prev_fpr)
        end
        prev_tpr = tpr
        prev_fpr = fpr
        idx += 1
    end
    auroc
end

"""
fast_auprc(scores::AbstractArray{Float64}, classes::AbstractArray{Symbol}, pos_label::Symbol)

Computes the area under the Precision-Recall curve using a lower
trapezoidal estimator, which is more accurate for skewed datasets.


K. Boyd, K. H. Eng, and C. D. Page, “Area under the Precision-Recall
Curve: Point Estimates and Confidence Intervals,” in Machine Learning
and Knowledge Discovery in Databases, H. Blockeel, K. Kersting,
S. Nijssen, and F. Železný, Eds. Springer Berlin Heidelberg, 2013,
pp. 451–466.
"""
function fast_auprc(scores::AbstractArray{Float64}, classes::AbstractArray{Symbol}, pos_label::Symbol)
    num_scores = length(scores)
    ordering = sortperm(scores, rev=true)
    labels = classes[ordering]

    num_pos, num_neg = 0, 0
    for label in labels
        if label == pos_label
            num_pos += 1
        else
            num_neg += 1
        end
    end

    tn, fn, tp, fp = 0, 0, num_pos, num_neg

    p = Array(Float64, num_scores)
    r = Array(Float64, num_scores)
    p[num_scores] = tp/(tp+fp)
    r[num_scores] = tp/(tp+fn)
    auprc, prev_r = 0.0, r[num_scores]
    pmin, pmax = p[num_scores], p[num_scores]

    # traverse scores from lowest to highest
    for i in num_scores-1:-1:1
        dtn = labels[i+1] != pos_label ? 1 : 0
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
