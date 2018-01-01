import pandas as pd
import os
def merge_bc_counts(data_dir, out_name):
    file_name = os.path.join(data_dir, "bc_count_{}.csv")
    df_matrix = None
    for i in range(1,5):
        df = pd.read_csv(file_name.format(i))
        df2 = df.loc[:,("gene", "barcodeid", "class", "counts_bin1", "counts_bin2")]
        df2.loc[:,"gene"] = [ "G%04d"%g for g in df2.loc[:,"gene"]]
        df2.loc[:,"barcodeid"] = [ "%s_%d"%(g,s) for g, s in zip(df2.loc[:,"gene"], df2.loc[:,"barcodeid"]) ]
        df2.loc[:,"counts_bin1"] -= 0.5
        df2.loc[:,"counts_bin2"] -= 0.5
        if df_matrix is None:
            df_matrix = df2.copy()
        else:
            df_matrix = df_matrix.merge(df2[["barcodeid","counts_bin1", "counts_bin2"]], on="barcodeid")
        LOW = "L{}".format(i)
        HIGH = "H{}".format(i)
        df_matrix.rename(columns={'counts_bin1': LOW, 'counts_bin2': HIGH}, inplace=True)

    df_matrix.rename(columns={"barcodeid":"sgRNA"}, inplace=True)
    cols = ["gene", "sgRNA", "class"] + ["L"+str(i) for i in range(1,5)] + ["H"+str(i) for i in range(1,5)]
    df_matrix[cols].to_csv(out_name, index=False)



def convert_to_mageck_format(in_name, out_name):
    df_org = pd.read_csv(in_name)
    cols = ["sgRNA", "gene"] + ["L"+str(i) for i in range(1,5)] + ["H"+str(i) for i in range(1,5)]
    df_org[cols].to_csv(out_name, index=False, sep="\t")
