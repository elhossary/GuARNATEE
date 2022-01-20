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
    parser.add_argument("--min_raw_height", default=10, type=float,
                        help="")
    parser.add_argument("--gff_out_dir", required=True, type=str,
                        help="")
    args = parser.parse_args()

    # load files
    gff = GFF(gff_paths=args.gffs)
    filtered_gff = gff.filter(anno_type="CDS", min_len=50)

    wig_info_df = pd.DataFrame([x.split(":") for x in args.wigs],
                               columns=["file_path", "strand", "condition", "replicate", "treatment"])
    wig_info_df["file_desc"] =\
        wig_info_df["strand"] + "_" + wig_info_df["condition"] + "_rep_" + wig_info_df["replicate"]
    wig_info_df.drop(["strand", "condition", "replicate"], inplace=True, axis=1)
    for file_disc in wig_info_df["file_desc"].unique():
        working_wigs = wig_info_df[wig_info_df["file_desc"] == file_disc].loc[:, ["file_path", "treatment"]]
        if working_wigs.shape[0] != 3:
            exit(1)
        working_pathes = dict(zip(working_wigs["treatment"], working_wigs["file_path"]))
        TEX_neg = Wiggle(os.path.abspath(working_pathes["TEX_neg"]))
        TEX_pos = Wiggle(os.path.abspath(working_pathes["TEX_pos"]))
        term = Wiggle(os.path.abspath(working_pathes["term"]))
        WindowSRNA(filtered_gff, TEX_neg, term)\
            .call_window_srna(args.min_len, args.max_len, args.read_length, args.threshold_factor, args.min_raw_height)
        WindowSRNA(filtered_gff, TEX_pos, term) \
            .call_window_srna(args.min_len, args.max_len, args.read_length, args.threshold_factor, args.min_raw_height)
    exit()


if __name__ == '__main__':

    main()