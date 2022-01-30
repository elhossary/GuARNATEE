import pandas as pd
import pybedtools as pybed
from io import StringIO
from winrna_lib.window_peaks import WindowPeaks
import os
from tqdm import tqdm


class WindowSRNA:

    def __init__(self,
                 five_end_wiggle,
                 three_end_wiggle):
        self.seqids = set.intersection(set(list(five_end_wiggle.signals.keys())),
                                       set(list(three_end_wiggle.signals.keys())))
        self.five_end_wiggle = five_end_wiggle.signals
        self.three_end_wiggle = three_end_wiggle.signals
        strands = set(list(five_end_wiggle.orientations.values()) +
                      list(three_end_wiggle.orientations.values()))
        self.strand = list(strands)[0] if len(strands) == 1 else None
        if self.strand is None:
            print("Error: Non-unified stranded wiggles passed, please unify the strand files")
            exit(1)
        self.srna_candidates = pd.DataFrame(columns=["seqid", "start", "end", "attributes"])
        self.gff_col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

    def call_window_srna(self, min_len: int, max_len: int, min_distance: int,
                         threshold_factor: float, min_height: float):
        for seqid in self.seqids:
            if seqid != "NC_002516.2":
                continue
            print(f"Processing wiggles of {seqid}")
            five_end_peaks_obj = WindowPeaks(self.five_end_wiggle[seqid],
                                             min_distance, threshold_factor, min_height,
                                             True if self.strand == "-" else False, "SS")
            three_end_peaks_obj = WindowPeaks(self.three_end_wiggle[seqid],
                                              min_distance, threshold_factor, min_height,
                                              True if self.strand == "+" else False, "TS")
            """
            five_end_peaks_obj.export_to_gff(
                "five.gff", seqid=seqid, strand=self.strand, anno_source="NA", anno_type="SS")
            three_end_peaks_obj.export_to_gff(
                "three.gff", seqid=seqid, strand=self.strand, anno_source="NA", anno_type="TS")
            """
            five_end_peaks_bed = pybed.BedTool(five_end_peaks_obj.get_bed_str(seqid), from_string=True).sort()
            three_end_peaks_bed = pybed.BedTool(three_end_peaks_obj.get_bed_str(seqid), from_string=True).sort()
            connected_peaks_df = \
                self.connect_sites(five_end_peaks_bed, three_end_peaks_bed, min_len, max_len) \
                    if self.strand == "+" else \
                    self.connect_sites(three_end_peaks_bed, five_end_peaks_bed, min_len, max_len)
            self.srna_candidates = pd.concat([self.srna_candidates, connected_peaks_df], ignore_index=True)

    @staticmethod
    def connect_sites(start_bed: pybed, end_bed: pybed, min_len: int, max_len: int) -> pd.DataFrame:
        base_columns = ["seqid", "start", "end", "attributes"]
        min_len -= 1
        max_len -= 1
        ret_df = pd.DataFrame(columns=base_columns)
        start_df = start_bed.to_dataframe(names=base_columns)
        end_df = end_bed.to_dataframe(names=base_columns)
        for row_id in tqdm(start_df.index, desc="==> Connecting 5' - 3' ends"):
            size_range = set(range(start_df.at[row_id, "start"] + min_len, start_df.at[row_id, "start"] + max_len, 1))
            tmp_df = end_df[(end_df["seqid"] == start_df.at[row_id, "seqid"]) &
                            (end_df["end"].isin(size_range))].copy()
            if tmp_df.empty:
                continue
            tmp_df["start"] = start_df.at[row_id, "start"]
            tmp_df["length"] = tmp_df["end"] - tmp_df["start"] + 1
            tmp_df["attributes"] = f'{start_df.at[row_id, "attributes"]};' + \
                                   tmp_df["attributes"] + ";length=" + tmp_df["length"].astype(str)
            tmp_df.drop(["length"], inplace=True, axis=1)
            ret_df = pd.concat([ret_df, tmp_df], ignore_index=True)
        ret_df.sort_values(["seqid", "start", "end"], inplace=True)
        ret_df.reset_index(inplace=True, drop=True)
        return ret_df

    def export_to_gff(self, out_path: str, anno_source="_", anno_type="_"):
        gff_df = self.srna_candidates.copy()
        non_gff_columns = [x for x in gff_df.columns.tolist() if x not in self.gff_col_names]
        gff_df["source"] = anno_source
        gff_df["type"] = anno_type
        gff_df["score"] = "."
        gff_df["strand"] = self.strand
        gff_df["phase"] = "."
        for i in gff_df.index:
            gff_df.at[i, "attributes"] = f'ID={anno_type}_{gff_df.at[i, "seqid"]}{self.strand}_{i}'\
                                         f';name={anno_type}_{gff_df.at[i, "seqid"]}{self.strand}_{i}'
            for col in non_gff_columns:
                gff_df.at[i, "attributes"] += f';{col}={gff_df.at[i, col]}'

        gff_df.drop(non_gff_columns, inplace=True, axis=1)
        gff_df = gff_df.reindex(columns=self.gff_col_names)
        gff_df.to_csv(os.path.abspath(out_path), index=False, sep="\t", header=False)
        print("GFF exported")
