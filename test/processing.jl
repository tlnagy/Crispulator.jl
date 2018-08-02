raw_data = OrderedDict{Symbol, DataFrame}()
low = DataFrame(
    class = [:negcontrol, :not, :not, :not],
    counts=[10, 100, 1000, 0],
    barcodeid=[1,2,3,4],
    gene=[1,2,2,3],
    behavior=[:linear,:linear,:linear,:linear]
)
high = DataFrame(
    class = [:not, :negcontrol, :not, :not],
    counts=[10, 100, 1000, 10],
    barcodeid=[2,1,3,4],
    gene=[2,1,2,3],
    behavior=[:linear,:linear,:linear,:linear]
)
raw_data[:low] = low
raw_data[:high] = high

guide_data, gene_data = differences_between_bins(raw_data)

total = sum([10, 100, 1000, 0] .+ 0.5)
@test guide_data[:freqs_low][1] == 10.5/total
@test guide_data[:freqs_low][4] == 0.5/total
@test guide_data[:rel_freqs_low][1] == 1.0 # only a single negative control should it should be divided by itself
@test guide_data[:rel_freqs_low][2] == 100.5/total/(10.5/total)

# make sure the log2 fold changes are correct
@test all(guide_data[:log2fc_high_div_low] .== log2.(guide_data[:rel_freqs_high]./guide_data[:rel_freqs_low]))

# > 2 bins
raw_data = OrderedDict{Symbol, DataFrame}()
raw_data[:low] = low

raw_data[:mid] = DataFrame(
    class = [:not, :negcontrol, :not, :not],
    counts=[10, 100, 1000, 100],
    barcodeid=[2,1,3,4],
    gene=[2,1,2,3],
    behavior=[:linear,:linear,:linear,:linear]
)

raw_data[:high] = high
guide_data, gene_data = differences_between_bins(raw_data)

# make sure all the log2fc comparisons are there
@test all(map(x->contains(==, names(guide_data), x), [:log2fc_mid_div_low, :log2fc_high_div_low, :log2fc_high_div_mid]))
