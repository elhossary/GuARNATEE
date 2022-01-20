import pandas as pd
import pybedtools as pybed
from io import StringIO
from winrna_lib.window_peaks import WindowPeaks
import os


class WindowSRNA:

    def __init__(self,
                 gff_df: pd.DataFrame,
                 five_end_wiggle,
                 three_end_wiggle):
        self.gff_df = gff_df
        self.seqids = set.intersection(set(list(five_end_wiggle.signals.keys())),
                                       set(list(three_end_wiggle.signals.keys())),
                                       set(gff_df["seqid"].unique().tolist()))
        self.five_end_wiggle = five_end_wiggle.signals
        self.three_end_wiggle = three_end_wiggle.signals
        strands = set(list(five_end_wiggle.orientations.values()) +
                      list(three_end_wiggle.orientations.values()))
        self.strand = list(strands)[0] if len(strands) == 1 else None
        if self.strand is None:
            print("Error: Non-unified stranded wiggles passed, please unify the strand files")
            exit(1)
        self.gff_df[self.gff_df["strand"] == self.strand].copy()
        self.gff_col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
        self.srna_candidates = pd.DataFrame(columns=["seqid", "start", "end", "length"])

    def call_window_srna(self, min_len: int, max_len: int, min_distance: int, threshold_factor=1.5, min_height=0.0):
        five_end_is_reversed = True if self.strand == "-" else False
        three_end_is_reversed = True if self.strand == "+" else False
        gff_temp = self.gff_df.copy()
        gff_temp["group_index"] = pd.Series(list(range(gff_temp.shape[0])))
        widows_df = gff_temp.loc[:, ["group_index", "start", "end", "attributes"]]
        windows = list(zip((widows_df["start"] - 1).tolist(), (widows_df["end"] - 1).tolist()))
        windows = dict(zip(gff_temp["group_index"].tolist(), windows))
        for seqid in self.seqids:
            print(f"Processing wiggles of {seqid}")
            five_end_peaks_obj = WindowPeaks(self.five_end_wiggle[seqid],
                                             windows,
                                             min_distance,
                                             threshold_factor,
                                             min_height,
                                             five_end_is_reversed)
            three_end_peaks_obj = WindowPeaks(self.three_end_wiggle[seqid],
                                              windows,
                                              min_distance,
                                              threshold_factor,
                                              min_height,
                                              three_end_is_reversed)
            five_end_peaks_bed = pybed.BedTool(five_end_peaks_obj.get_bed_str(seqid), from_string=True).sort()
            three_end_peaks_bed = pybed.BedTool(three_end_peaks_obj.get_bed_str(seqid), from_string=True).sort()
            if self.strand == "+":
                connected_peaks = five_end_peaks_bed.closest(three_end_peaks_bed, d=True, D="ref", io=True, iu=True)
            else:
                connected_peaks = three_end_peaks_bed.closest(five_end_peaks_bed, d=True, D="ref", io=True, iu=True)
            connected_peaks_df = pd.read_csv(StringIO(str(connected_peaks)),
                                             names=["seqid", "start", "drop1", "drop2", "end", "drop3", "drop4"],
                                             sep="\t")
            connected_peaks_df["length"] = connected_peaks_df["end"] - connected_peaks_df["start"] + 1
            connected_peaks_df = connected_peaks_df[connected_peaks_df["drop4"].isin(range(min_len, max_len - 1, 1))]
            print(connected_peaks_df.to_string())
            drop_columns = [x for x in connected_peaks_df.columns.tolist() if "drop" in x]
            connected_peaks_df.drop(drop_columns, axis=1, inplace=True)
            connected_peaks_df["length"] = connected_peaks_df["end"] - connected_peaks_df["start"] + 1
            #connected_peaks_df = connected_peaks_df[connected_peaks_df["length"].between(min_len, max_len)]
            connected_peaks_df.drop_duplicates(inplace=True)
            self.srna_candidates = self.srna_candidates.append(connected_peaks_df, ignore_index=True)

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
