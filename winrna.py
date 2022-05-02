import sys
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
import pybedtools as pybed
import configparser

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
        "--config_file", default=f"{os.path.dirname(__file__)}/config.cfg", type=str, help="Configuration file"
    )
    parser.add_argument("--out_dir", required=True, type=str, help="")
    args = parser.parse_args()

    # load files
    conf_dict = {}
    conf_file_path = os.path.abspath(args.config_file)
    conf_parse = configparser.ConfigParser(strict=True)
    print(f"Loading config file: {conf_file_path}")
    conf_parse.read(conf_file_path)
    for sect in conf_parse.sections():
        conf_dict.update(dict(conf_parse.items(sect)), )

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
    stats_df = pd.DataFrame()
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
                sys.exit(1)
            working_pathes = dict(
                zip(working_wigs["treatment"], working_wigs["file_path"])
            )
            treated_srnas_df, treated_stats_df = _call_srnas(
                working_pathes["TEX_pos"], working_pathes["term"], conf_dict
            )
            treated_stats_df["file_desc"] = desc
            treated_stats_df["TSS_lib_type"] = "treated"
            treated_stats_df["strand"] = strand

            treated_srnas_df["strand"] = strand_sign
            control_srnas_df, control_stats_df = _call_srnas(
                working_pathes["TEX_neg"], working_pathes["term"], conf_dict
            )
            control_stats_df["file_desc"] = desc
            control_stats_df["TSS_lib_type"] = "control"
            control_stats_df["strand"] = strand

            control_srnas_df["strand"] = strand_sign
            tmp_df1 = pd.concat([tmp_df1, treated_srnas_df], ignore_index=True)
            tmp_df2 = pd.concat([tmp_df2, control_srnas_df], ignore_index=True)
            stats_df = pd.concat([stats_df, treated_stats_df, control_stats_df], ignore_index=True)

        tmp_df1 = Helpers.get_gff_df(
            tmp_df1, anno_source="WinRNA", anno_type="candidate", new_id=True
        )
        tmp_df2 = Helpers.get_gff_df(
            tmp_df2, anno_source="WinRNA", anno_type="candidate", new_id=True
        )
        tmp_df1["condition"] = desc
        tmp_df2["condition"] = desc
        tmp_df1 = Helpers.warp_non_gff_columns(RNAClassifier(gff_obj, tmp_df1, fastas, conf_dict).classes)
        tmp_df2 = Helpers.warp_non_gff_columns(RNAClassifier(gff_obj, tmp_df2, fastas, conf_dict).classes)
        tmp_df1, tmp_df2 = DifferentialClassifier(
            {"TEX_pos": tmp_df1, "TEX_neg": tmp_df2}
        ).score_similarity()
        tmp_df1["replicate"] = desc
        tmp_df2["replicate"] = desc
        export_df = pd.concat([export_df, tmp_df1, tmp_df2], ignore_index=True)
    export_df = Helpers.warp_non_gff_columns(export_df)
    export_pb = pybed.BedTool.from_dataframe(export_df).sort().cluster(s=True, d=0)
    export_df = export_pb.to_dataframe(names=gff_obj.column_names + ["cluster_id"])
    clusters = export_df["cluster_id"].unique().size
    export_df = Helpers.warp_non_gff_columns(export_df)
    export_df.sort_values(["seqid", "start", "end"], inplace=True)
    export_df.reset_index(inplace=True, drop=True)
    export_df["source"] = "WinRNA"
    print(f"Total {export_df.shape[0]} candidates in {clusters} unique regions to be exported")
    # Exports
    # ==> GFFs
    for seqid_group, seqids in seqid_groups.items():
        export_df[export_df["seqid"].isin(seqids)].to_csv(
            os.path.abspath(f"{args.out_dir}/{seqid_group}_candidates.gff"),
            index=False,
            sep="\t",
            header=False)
    # Stats
    for i in stats_df.index:
        for k, v in seqid_groups.items():
            if stats_df.at[i, "seqid"] in v:
                stats_df.at[i, "Organism"] = k
    stats_df = stats_df.groupby(["Organism", "file_desc", "TSS_lib_type"], as_index=False).agg({"TSS_lib_windows_count": "sum",
                                                                                                "TSS_lib_peaks_count": "sum",
                                                                                                "TTS_lib_windows_count": "sum",
                                                                                                "TTS_lib_peaks_count": "sum",
                                                                                                "peaks_connections_count": "sum"})
    stats_df.to_csv(os.path.abspath(f"{args.out_dir}/stats.tsv"), sep='\t', index=False)

    # ==> Excel
    drop_cols = ["source", "type", "score", "phase", "id",
                 "ss_id", "ss_diff_height", "ss_height", "ss_upstream", "ss_downstream",
                 "ts_id", "ts_diff_height", "ts_height", "ts_upstream", "ts_downstream"]

    export_df = Helpers.expand_attributes_to_columns(export_df)
    export_df.drop(columns=drop_cols, inplace=True)
    rename_cols = {c: c.replace("_", " ") for c in export_df.columns}

    classes_groups = {"intergenic": [], "ORF_int": [], "others": []}
    for cls in export_df["annotation_class"].unique():
        if ("cross" in cls and "ncRNA" in cls) or "antisense_to" in cls:
            classes_groups["others"] += [cls]
        elif "ncRNA" in cls or "intergenic" in cls:
            classes_groups["intergenic"] += [cls]
        elif "ORF_int" in cls:
            classes_groups["ORF_int"] += [cls]
        else:
            classes_groups["others"] += [cls]

    ### ORF int Stats
    orf_int_stats_df = export_df[export_df["sub_class"] != ""].copy()
    orf_int_stats_df["Specie"] = ""
    for seqid_group, seqids in seqid_groups.items():
        orf_int_stats_df.loc[orf_int_stats_df["seqid"].isin(seqids), "Specie"] = seqid_group
    orf_int_stats_df = orf_int_stats_df.value_counts(subset=["Specie", "lib_type", "replicate", "sub_class"]).reset_index()
    orf_int_stats_df.columns = ["Specie", "lib_type", "replicate", "sub_class", 'count']
    with pd.ExcelWriter(os.path.abspath(f"{args.out_dir}/ORF_int_stats.xlsx"), engine="openpyxl") as writer:
        for sp in orf_int_stats_df["Specie"].unique():
            tmp_stat_df = orf_int_stats_df[orf_int_stats_df["Specie"] == sp].copy()
            #tmp_stat_df.drop(columns=[sp], inplace=True)
            tmp_stat_df.to_excel(excel_writer=writer, sheet_name=sp, index=False, header=True, na_rep="", verbose=True)

    ### Tables
    for seqid_group, seqids in seqid_groups.items():
        with pd.ExcelWriter(os.path.abspath(f"{args.out_dir}/{seqid_group}_candidates.xlsx"), engine="openpyxl") \
                as writer:
            for classes_group, classes in classes_groups.items():
                seqid_type_df = export_df[(export_df["seqid"].isin(seqids)) &
                                          (export_df["annotation_class"].isin(classes))].copy()
                if seqid_type_df.empty:
                    continue
                seqid_type_df.reset_index(inplace=True, drop=True)
                seqid_type_df.replace("", np.nan, inplace=True)
                seqid_type_df.dropna(how='all', axis=1, inplace=True)
                seqid_type_df.rename(columns=rename_cols, inplace=True)
                seqid_type_df.to_excel(
                    excel_writer=writer,
                    sheet_name=f"{classes_group}",
                    index=True,
                    header=True,
                    na_rep="",
                    verbose=True,
                    index_label="index")
    sys.exit(0)


def _call_srnas(five_end_path, three_end_path, conf_dict):
    srnas = WindowSRNA(Wiggle(five_end_path), Wiggle(three_end_path))
    srnas.call_window_srna(conf_dict)
    return srnas.srna_candidates, srnas.log_df


if __name__ == "__main__":
    main()
