function compute_roc(df::DataFrame, score_column::Symbol, numsteps::Int64)
    tprs = zeros(numsteps)
    fprs = zeros(numsteps)
    auroc = compute_roc!(Vector(df[score_column]), Vector(df[:class]) , numsteps, tprs, fprs)
    return auroc, tprs, fprs
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
        tp, fp, tn, fn = 0,0,0,0

        for i in 1:num_scores
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
        tp, fp, tn, fn = 0,0,0,0

        for i in 1:num_scores
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
