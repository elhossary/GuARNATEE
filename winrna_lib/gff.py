import glob
import os.path
import pandas as pd


class GFF:

    def __init__(self, gff_paths: list):
        self.gff_paths = gff_paths
        self.column_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
        self.gff_df = pd.DataFrame(columns=self.column_names)
        self.parse()

    def parse(self):
        print("=> Parsing input GFF file")
        parsed_paths = []
        for item in self.gff_paths:
            for sub_item in glob.glob(item):
                parsed_paths.append(os.path.abspath(sub_item))
        for gff_path in parsed_paths:
            self.gff_df = self.gff_df.append(
                pd.read_csv(gff_path, names=self.column_names, comment="#", sep="\t"),
                ignore_index=True)
        self.gff_df.reset_index(drop=True, inplace=True)

    def filter(self, anno_type=None, min_len=0, max_len=0, inplace=False) -> pd.DataFrame:
        gff_df = self.gff_df.copy()
        if anno_type is not None:
            types_list = str(anno_type).split(" ")
            gff_df = gff_df[gff_df["type"].isin(types_list)]
        if min_len > 0 or max_len > 0:
            gff_df["length"] = gff_df["end"] + gff_df["start"] - 1
            if min_len > 0:
                gff_df = gff_df[gff_df["length"] >= max_len]
            if max_len > 0:
                gff_df = gff_df[gff_df["length"] <= max_len]
            gff_df.drop(["length"], inplace=True, axis=1)
        if inplace:
            self.gff_df = gff_df
            return self.gff_df
        else:
            return gff_df

    @staticmethod
    def write_to_gff(peaks_df: pd.DataFrame, out_path: str) -> None:
        strand_func = lambda x: "F" if x == "+" else "R"
        wrap_attr_func = lambda x: ";".join([f"{k}={v}" for k, v in x.items()])
        col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
        peaks_df["source"] = ""
        peaks_df["type"] = ""
        peaks_df["score"] = "."
        peaks_df["phase"] = "."
        peaks_df = peaks_df.round(2)
        atrr_cols = [col for col in peaks_df.columns.tolist() if col not in col_names]
        peaks_df["attributes"] = ""
        for i in peaks_df.index:
            peaks_df.at[i, "attributes"] = \
                f"ID={peaks_df.at[i, 'seqid']}{strand_func(peaks_df.at[i, 'strand'])}_{i};" + \
                f"name={peaks_df.at[i, 'seqid']}{strand_func(peaks_df.at[i, 'strand'])}_{i};" + \
                f"{wrap_attr_func(peaks_df.loc[i, atrr_cols])}"

        peaks_df.drop(atrr_cols, inplace=True, axis=1)
        peaks_df = peaks_df.reindex(col_names, axis="columns")
        out_gff = f"{wiggle_path}.gff"
        if args.output_dir is not None:
            wiggle_basename = os.path.basename(wiggle_path)
            out_gff = f"{os.path.dirname(args.output_dir)}/{wiggle_basename}.gff"
        peaks_df.to_csv(out_gff, sep="\t", header=False, index=False)
