from winrna_lib.wiggle import Wiggle
from winrna_lib.window_srna import WindowSRNA
from winrna_lib.gff import GFF
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
    parser.add_argument("--threshold_factor", default=1.5, type=float,
                        help="")
    parser.add_argument("--min_raw_height", default=0, type=float,
                        help="")
    parser.add_argument("--gff_out_dir", required=True, type=str,
                        help="")
    args = parser.parse_args()

    # load files
    gff = GFF(gff_paths=args.gffs)
    filtered_gff = gff.filter(anno_type="CDS", min_len=50)
    wig_info_df = pd.DataFrame([x.split(":") for x in args.wigs],
                               columns=["file_path", "strand", "condition", "replicate", "treatment"])
    wig_info_df["file_desc"] = wig_info_df["condition"] + "_rep_" + wig_info_df["replicate"]

    for desc in wig_info_df["file_desc"].unique().tolist():
        tmp_df1 = pd.DataFrame()
        tmp_df2 = pd.DataFrame()
        for strand in ["f", "r"]:
            strand_sign = "+" if strand == "f" else "-"
            stranded_filtered_gff = filtered_gff[filtered_gff["strand"] == strand_sign]
            working_wigs = wig_info_df[(wig_info_df["strand"] == strand) &
                                       (wig_info_df["file_desc"] == desc)].loc[:, ["file_path", "treatment"]]
            if working_wigs.shape[0] != 3:
                exit(1)
            working_pathes = dict(zip(working_wigs["treatment"], working_wigs["file_path"]))
            treated_srnas_df = \
                _call_srnas(working_pathes["TEX_pos"], working_pathes["term"], stranded_filtered_gff, args)
            treated_srnas_df["strand"] = strand_sign
            control_srnas_df = \
                _call_srnas(working_pathes["TEX_neg"], working_pathes["term"], stranded_filtered_gff, args)
            control_srnas_df["strand"] = strand_sign
            tmp_df1 = tmp_df1.append(treated_srnas_df, ignore_index=True)
            tmp_df2 = tmp_df2.append(control_srnas_df, ignore_index=True)

        export_to_gff(tmp_df1,
                      f"{os.path.dirname(args.gff_out_dir)}/TEX_pos_{desc}.gff", "WinRNA", "ORF_int")
        export_to_gff(tmp_df2,
                      f"{os.path.dirname(args.gff_out_dir)}/TEX_neg_{desc}.gff", "WinRNA", "ORF_int")
    exit(0)


def _call_srnas(five_end_path, three_end_path, stranded_filtered_gff, args):
    srnas = WindowSRNA(stranded_filtered_gff, Wiggle(five_end_path), Wiggle(three_end_path))
    srnas.call_window_srna(
        args.min_len, args.max_len, args.read_length, args.threshold_factor, args.min_raw_height)
    return srnas.srna_candidates


def export_to_gff(gff_df: pd.DataFrame, out_path: str, anno_source="_", anno_type="_"):
    gff_col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    non_gff_columns = [x for x in gff_df.columns.tolist() if x not in gff_col_names]
    gff_df["source"] = anno_source
    gff_df["type"] = anno_type
    gff_df["score"] = "."
    gff_df["phase"] = "."
    for i in gff_df.index:
        strand = "F" if gff_df.at[i, "strand"] == "+" else "R"
        gff_df.at[i, "attributes"] = f'ID={anno_type}_{gff_df.at[i, "seqid"]}{strand}_{i}'\
                                     f';name={anno_type}_{gff_df.at[i, "seqid"]}{strand}_{i}' \
                                     f';{gff_df.at[i, "attributes"]}'

    gff_df.drop(non_gff_columns, inplace=True, axis=1)
    gff_df = gff_df.reindex(columns=gff_col_names)
    gff_df.sort_values(["seqid", "start", "end"], inplace=True)
    gff_df.to_csv(os.path.abspath(out_path), index=False, sep="\t", header=False)
    print("GFF exported")

def export_to_table():
    pass

if __name__ == '__main__':

    main()
