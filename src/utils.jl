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
compute_auprc(scores::Vector{Float64}, classes::Vector{Symbol})

Computes the area under the Precision-Recall curve using a lower
trapezoidal estimator, which is more accurate for skewed datasets.


K. Boyd, K. H. Eng, and C. D. Page, “Area under the Precision-Recall
Curve: Point Estimates and Confidence Intervals,” in Machine Learning
and Knowledge Discovery in Databases, H. Blockeel, K. Kersting,
S. Nijssen, and F. Železný, Eds. Springer Berlin Heidelberg, 2013,
pp. 451–466.
"""
function compute_auprc(scores::Vector{Float64}, classes::Vector{Symbol})

    num_scores = length(scores)
    thresholds = sort(scores, rev=true)
    auprc, prev_recall, pmin,pmax = 0.0, -1.0, -Inf, Inf
    precision = Array(Float64, num_scores)
    recall = Array(Float64, num_scores)
    tp, fp, tn, fn = _binary_classifier(scores, classes, thresholds[num_scores])

    precision[num_scores] = tp/(tp+fp)
    recall[num_scores] = tp/(tp+fn)
    prev_recall = recall[num_scores]
    pmin = precision[num_scores]

    @inbounds for i in num_scores-1:-1:1

        tp, fp, tn, fn = _binary_classifier(scores, classes, thresholds[i])

        precision[i] = tp/(tp+fp)
        recall[i] = tp/(tp+fn)

        if recall[i] == prev_recall
            pmax = precision[i]
        else
            pmin = precision[i]
            auprc += (pmin + pmax)/2*(prev_recall - recall[i])
            prev_recall = recall[i]
        end
    end
    auprc, precision, recall
end
