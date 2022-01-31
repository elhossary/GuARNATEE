from winrna_lib.wiggle import Wiggle
from winrna_lib.window_srna import WindowSRNA
from winrna_lib.gff import GFF
from winrna_lib.rna_classifier import RNAClassifier
import argparse
import pandas as pd
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gffs", required=True, type=str, nargs="+",
                        help="GFF files (space separated)")
    parser.add_argument("--wigs", required=True, type=str, nargs="+",
                        help="Wiggle files (space separated)")
    parser.add_argument("--min_len", default=50, type=int,
                        help="Minimum allowed annotation length")
    parser.add_argument("--max_len", default=300, type=int,
                        help="Maximum allowed annotation length")
    parser.add_argument("--read_length", default=75, type=int,
                        help="")
    parser.add_argument("--threshold_factor", default=1, type=float,
                        help="")
    parser.add_argument("--min_raw_height", default=10, type=float,
                        help="")
    parser.add_argument("--min_bg_fold_change", default=1, type=float,
                        help="")
    parser.add_argument("--gff_out_dir", required=True, type=str,
                        help="")
    args = parser.parse_args()

    # load files
    gff_df = GFF(gff_paths=args.gffs).gff_df
    wig_info_df = pd.DataFrame([x.split(":") for x in args.wigs],
                               columns=["file_path", "strand", "condition", "replicate", "treatment"])
    wig_info_df["file_desc"] = wig_info_df["condition"] + "_rep_" + wig_info_df["replicate"]

    for desc in wig_info_df["file_desc"].unique().tolist():
        tmp_df1 = pd.DataFrame()
        tmp_df2 = pd.DataFrame()
        for strand in ["f", "r"]:
            strand_sign = "+" if strand == "f" else "-"
            working_wigs = wig_info_df[(wig_info_df["strand"] == strand) &
                                       (wig_info_df["file_desc"] == desc)].loc[:, ["file_path", "treatment"]]
            if working_wigs.shape[0] not in [3]:
                exit(1)
            working_pathes = dict(zip(working_wigs["treatment"], working_wigs["file_path"]))
            treated_srnas_df = \
                _call_srnas(working_pathes["TEX_pos"], working_pathes["term"], args)
            treated_srnas_df["strand"] = strand_sign
            control_srnas_df = \
                _call_srnas(working_pathes["TEX_neg"], working_pathes["term"], args)
            control_srnas_df["strand"] = strand_sign
            tmp_df1 = pd.concat([tmp_df1, treated_srnas_df], ignore_index=True)
            tmp_df2 = pd.concat([tmp_df2, control_srnas_df], ignore_index=True)
        tmp_df1 = to_gff_df(tmp_df1, "WinRNA")
        tmp_df2 = to_gff_df(tmp_df2, "WinRNA")
        tmp_df1 = RNAClassifier(gff_df, tmp_df1).classes
        tmp_df2 = RNAClassifier(gff_df, tmp_df2).classes

        # Exports
        tmp_df1.to_csv(os.path.abspath(f"{os.path.dirname(args.gff_out_dir)}/TEX_pos_{desc}.gff"),
                       index=False, sep="\t", header=False)
        tmp_df2.to_csv(os.path.abspath(f"{os.path.dirname(args.gff_out_dir)}/TEX_neg_{desc}.gff"),
                       index=False, sep="\t", header=False)
        to_table_df(tmp_df1) \
            .to_csv(os.path.abspath(f"{os.path.dirname(args.gff_out_dir)}/TEX_pos_{desc}.tsv"),
                    index=True, sep="\t", header=True)
        to_table_df(tmp_df2) \
            .to_csv(os.path.abspath(f"{os.path.dirname(args.gff_out_dir)}/TEX_neg_{desc}.tsv"),
                    index=True, sep="\t", header=True)
    exit(0)


def _call_srnas(five_end_path, three_end_path, args):
    srnas = WindowSRNA(Wiggle(five_end_path), Wiggle(three_end_path))
    srnas.call_window_srna(
        args.min_len, args.max_len, args.read_length, args.threshold_factor, args.min_raw_height)
    return srnas.srna_candidates


def to_gff_df(gff_df: pd.DataFrame, anno_source="_") -> pd.DataFrame:
    anno_type = "candidate"
    gff_col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    non_gff_columns = [x for x in gff_df.columns.tolist() if x not in gff_col_names]
    gff_df["source"] = anno_source
    gff_df["type"] = anno_type
    gff_df["score"] = "."
    gff_df["phase"] = "."
    for i in gff_df.index:
        strand = "F" if gff_df.at[i, "strand"] == "+" else "R"
        gff_df.at[i, "attributes"] = f'ID={gff_df.at[i, "seqid"]}{strand}_{anno_type}_{i}'\
                                     f';name={gff_df.at[i, "seqid"]}{strand}_{anno_type}_{i}' \
                                     f';{gff_df.at[i, "attributes"]}'
        for col in non_gff_columns:
            gff_df.at[i, "attributes"] += f";{col}={gff_df.at[i, col]}"
    gff_df.drop(non_gff_columns, inplace=True, axis=1)
    gff_df = gff_df.reindex(columns=gff_col_names)
    gff_df.sort_values(["seqid", "start", "end"], inplace=True)
    return gff_df


def to_table_df(df):
    df = expand_attributes_to_columns(df)
    df.sort_values(["seqid", "start", "end"], inplace=True)
    df.reset_index(inplace=True, drop=True)
    return df


def expand_attributes_to_columns(df):
    if "attributes" in df.columns:
        for i in df.index:
            attributes = df.at[i, "attributes"].split(";")
            for attr in attributes:
                k, v = attr.split("=")
                df.at[i, k] = v
        df.drop(["attributes"], inplace=True, axis=1)
    return df



if __name__ == '__main__':

    main()
