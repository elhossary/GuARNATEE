import sys
import os
import pandas as pd
import pybedtools as pybed
from tqdm import tqdm
from winrna_lib.window_peaks import WindowPeaks
from winrna_lib.helpers import Helpers


class WindowSRNA:
    """
    This class is build upon WindowPeaks class,
    which combines peaks called from WindowPeaks into WindowSRNA to form annotation candidates
    """

    def __init__(self, five_end_wiggle, three_end_wiggle):
        self.seqids = set.intersection(
            set(list(five_end_wiggle.signals.keys())),
            set(list(three_end_wiggle.signals.keys())),
        )
        self.five_end_wiggle = five_end_wiggle.signals
        self.three_end_wiggle = three_end_wiggle.signals
        strands = set(
            list(five_end_wiggle.orientations.values())
            + list(three_end_wiggle.orientations.values())
        )
        self.strand = list(strands)[0] if len(strands) == 1 else None
        if self.strand is None:
            print(
                "Error: Non-unified stranded wiggles passed, please unify the strand files"
            )
            sys.exit(1)
        self.srna_candidates = pd.DataFrame(
            columns=["seqid", "start", "end", "attributes"]
        )
        self.gff_col_names = [
            "seqid",
            "source",
            "type",
            "start",
            "end",
            "score",
            "strand",
            "phase",
            "attributes",
        ]
        self.log_df = pd.DataFrame(
            columns=["seqid",
                     "TSS_lib_windows_count", "TSS_lib_peaks_count",
                     "TTS_lib_windows_count", "TTS_lib_peaks_count",
                     "peaks_connections_count"])

    def call_window_srna(
        self, conf_dict: dict, thres_factor) -> None:
        for seqid in self.seqids:
            """
            # for testing purpose only
            if seqid != "NC_002516.2":
                continue
            """
            print(f"=> Calling 5' ends for SeqID: {seqid}")
            five_end_peaks_obj = WindowPeaks(
                self.five_end_wiggle[seqid],
                conf_dict["min_distance"],
                conf_dict["min_height"],
                conf_dict["min_step_factor"],
                bool(self.strand == "-"),
                "SS",
                thres_factor
            )
            print(f"=> Calling 3' ends for SeqID: {seqid}")
            three_end_peaks_obj = WindowPeaks(
                self.three_end_wiggle[seqid],
                conf_dict["min_distance"],
                conf_dict["min_height"],
                conf_dict["min_step_factor"],
                bool(self.strand == "+"),
                "TS",
                thres_factor
            )

            five_end_peaks_str = five_end_peaks_obj.get_bed_str(seqid)
            three_end_peaks_str = three_end_peaks_obj.get_bed_str(seqid)
            if five_end_peaks_str is None or three_end_peaks_str is None:
                continue
            five_end_peaks_bed = pybed.BedTool(
                five_end_peaks_str, from_string=True
            ).sort()
            three_end_peaks_bed = pybed.BedTool(
                three_end_peaks_str, from_string=True
            ).sort()
            connected_peaks_df = (
                self.connect_sites(
                    five_end_peaks_bed, three_end_peaks_bed, int(conf_dict["min_len"]), int(conf_dict["max_len"])
                )
                if self.strand == "+"
                else self.connect_sites(
                    three_end_peaks_bed, five_end_peaks_bed, int(conf_dict["min_len"]), int(conf_dict["max_len"])
                )
            )
            tmp_dict = {"seqid": seqid,
                        "TSS_lib_windows_count": len(five_end_peaks_obj.windows),
                        "TSS_lib_peaks_count": five_end_peaks_obj.peaks_arr.size,
                        "TTS_lib_windows_count": len(three_end_peaks_obj.windows),
                        "TTS_lib_peaks_count": three_end_peaks_obj.peaks_arr.size,
                        "peaks_connections_count": connected_peaks_df.shape[0],
                        "TSS_peaks_thresholds": five_end_peaks_obj.peaks_arr[:, 1].tolist(),
                        "TTS_peaks_thresholds": three_end_peaks_obj.peaks_arr[:, 1].tolist()}
            self.log_df = pd.concat([self.log_df, pd.DataFrame([tmp_dict], columns=tmp_dict.keys())], ignore_index=True)
            self.srna_candidates = pd.concat(
                [self.srna_candidates, connected_peaks_df], ignore_index=True
            )

    @staticmethod
    def connect_sites(
        start_bed: pybed, end_bed: pybed, min_len: int, max_len: int
    ) -> pd.DataFrame:
        base_columns = ["seqid", "start", "end", "attributes"]
        min_len -= 1
        max_len -= 1
        ret_df = pd.DataFrame(columns=base_columns)
        start_df = start_bed.to_dataframe(names=base_columns)
        end_df = end_bed.to_dataframe(names=base_columns)
        for row_id in tqdm(start_df.index, desc="==> Connecting 5' - 3' ends", bar_format='{desc} |{bar:20}| {percentage:3.0f}%'):
            size_range = set(
                range(
                    start_df.at[row_id, "start"] + min_len,
                    start_df.at[row_id, "start"] + max_len,
                    1,
                )
            )
            tmp_df = end_df[
                (end_df["seqid"] == start_df.at[row_id, "seqid"])
                & (end_df["end"].isin(size_range))
            ].copy()
            if tmp_df.empty:
                continue
            tmp_df["start"] = start_df.at[row_id, "start"]
            tmp_df["length"] = tmp_df["end"] - tmp_df["start"] + 1
            tmp_df["attributes"] = (
                f'{start_df.at[row_id, "attributes"]};'
                + tmp_df["attributes"]
                + ";length="
                + tmp_df["length"].astype(str)
            )
            tmp_df.drop(["length"], inplace=True, axis=1)
            ret_df = pd.concat([ret_df, tmp_df], ignore_index=True)
        ret_df.sort_values(["seqid", "start", "end"], inplace=True)
        ret_df.reset_index(inplace=True, drop=True)
        return ret_df

    def export_to_gff(self, out_path: str, anno_source="NA", anno_type="NA"):
        gff_df = Helpers.get_gff_df(
            self.srna_candidates, anno_source=anno_source, anno_type=anno_type
        )
        gff_df.to_csv(os.path.abspath(out_path), index=False, sep="\t", header=False)
        print("GFF exported")
