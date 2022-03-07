import numpy as np

from winrna_lib.window_srna import WindowSRNA
from winrna_lib.rna_classifier import RNAClassifier
from winrna_lib.differential_classifier import DifferentialClassifier
from winrna_lib.helpers import Helpers
from winrna_lib.wiggle import Wiggle
from winrna_lib.gff import GFF
import argparse
import pandas as pd
import os
from winrna_lib.fasta import Fasta


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gffs", required=True, type=str, nargs="+", help="GFF files (space separated)"
    )
    parser.add_argument(
        "--fastas",
        required=True,
        type=str,
        nargs="+",
        help="Fasta files (space separated)",
    )
    parser.add_argument(
        "--wigs",
        required=True,
        type=str,
        nargs="+",
        help="Wiggle files (space separated)",
    )
    parser.add_argument(
        "--min_len", default=50, type=int, help="Minimum allowed annotation length"
    )
    parser.add_argument(
        "--max_len", default=350, type=int, help="Maximum allowed annotation length"
    )
    parser.add_argument("--read_length", default=75, type=int, help="")
    parser.add_argument("--min_raw_height", default=10, type=float, help="")
    parser.add_argument("--out_dir", required=True, type=str, help="")
    args = parser.parse_args()

    # load files
    gff_obj = GFF(gff_paths=args.gffs)
    fastas = Fasta(fasta_paths=args.fastas)
    seqid_groups = fastas.organism_seqid_groups
    wig_info_df = pd.DataFrame(
        [x.split(":") for x in args.wigs],
        columns=["file_path", "strand", "condition", "replicate", "treatment"],
    )
    wig_info_df["file_desc"] = (
        wig_info_df["condition"] + "_rep_" + wig_info_df["replicate"]
    )
    export_df = pd.DataFrame()
    for desc in wig_info_df["file_desc"].unique().tolist():
        tmp_df1 = pd.DataFrame()
        tmp_df2 = pd.DataFrame()
        for strand in ["f", "r"]:
            strand_sign = "+" if strand == "f" else "-"
            working_wigs = wig_info_df[
                (wig_info_df["strand"] == strand) & (wig_info_df["file_desc"] == desc)
            ].loc[:, ["file_path", "treatment"]]
            if working_wigs.shape[0] not in [3]:
                print("Error: non-uniformed wiggles passed")
                exit(1)
            working_pathes = dict(
                zip(working_wigs["treatment"], working_wigs["file_path"])
            )
            treated_srnas_df = _call_srnas(
                working_pathes["TEX_pos"], working_pathes["term"], args
            )
            treated_srnas_df["strand"] = strand_sign
            control_srnas_df = _call_srnas(
                working_pathes["TEX_neg"], working_pathes["term"], args
            )
            control_srnas_df["strand"] = strand_sign
            tmp_df1 = pd.concat([tmp_df1, treated_srnas_df], ignore_index=True)
            tmp_df2 = pd.concat([tmp_df2, control_srnas_df], ignore_index=True)
        tmp_df1 = Helpers.get_gff_df(
            tmp_df1, anno_source="WinRNA", anno_type="candidate", new_id=True
        )
        tmp_df2 = Helpers.get_gff_df(
            tmp_df2, anno_source="WinRNA", anno_type="candidate", new_id=True
        )
        tmp_df1 = Helpers.warp_non_gff_columns(RNAClassifier(gff_obj, tmp_df1, fastas).classes)
        tmp_df2 = Helpers.warp_non_gff_columns(RNAClassifier(gff_obj, tmp_df2, fastas).classes)
        tmp_df1, tmp_df2 = DifferentialClassifier(
            {"TEX_pos": tmp_df1, "TEX_neg": tmp_df2}
        ).score_similarity()
        tmp_df1["replicate"] = desc
        tmp_df2["replicate"] = desc
        export_df = pd.concat([export_df, tmp_df1, tmp_df2], ignore_index=True)
    export_df.sort_values(["seqid", "start", "end"], inplace=True)
    export_df.reset_index(inplace=True, drop=True)
    export_df["source"] = "WinRNA"

    # Exports
    ## GFFs
    for seqid_group, seqids in seqid_groups.items():
        Helpers.warp_non_gff_columns(export_df[export_df["seqid"].isin(seqids)]).to_csv(
            os.path.abspath(f"{args.out_dir}/{seqid_group}_candidates.gff"),
            index=False,
            sep="\t",
            header=False)
    export_df = Helpers.expand_attributes_to_columns(export_df)
    rename_cols = {c: c.replace("_", " ") for c in export_df.columns}
    for anno_type in export_df["annotation_class"].unique():
        type_df = export_df[export_df["annotation_class"] == anno_type]
        with pd.ExcelWriter(os.path.abspath(f"{args.out_dir}/{anno_type}_candidates.xlsx"), engine="openpyxl") as writer:
            for seqid_group, seqids in seqid_groups.items():
                export_tmp_df = type_df[type_df["seqid"].isin(seqids)].copy()
                export_tmp_df.reset_index(inplace=True, drop=True)
                export_tmp_df.replace("", np.nan, inplace=True)
                export_tmp_df.dropna(how='all', axis=1, inplace=True)
                export_tmp_df.rename(columns=rename_cols, inplace=True)
                export_tmp_df.to_excel(
                    excel_writer=writer,
                    sheet_name=f"{seqid_group}_{anno_type}",
                    index=True,
                    header=True,
                    na_rep="",
                    verbose=True)

    exit(0)


def _call_srnas(five_end_path, three_end_path, args):
    srnas = WindowSRNA(Wiggle(five_end_path), Wiggle(three_end_path))
    srnas.call_window_srna(
        args.min_len,
        args.max_len,
        args.read_length,
        args.min_raw_height,
    )
    return srnas.srna_candidates


if __name__ == "__main__":
    main()
