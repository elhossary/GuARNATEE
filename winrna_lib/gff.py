import glob
import os.path
import pandas as pd


class GFF:

    def __init__(self, gff_paths: list):
        self.gff_paths = gff_paths
        self.column_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
        self.gff_df = pd.DataFrame(columns=self.column_names)
        self.regions = pd.DataFrame()
        self.parse()

    def parse(self):
        print("=> Parsing input GFF file")
        parsed_paths = []
        for item in self.gff_paths:
            for sub_item in glob.glob(item):
                parsed_paths.append(os.path.abspath(sub_item))
        for gff_path in parsed_paths:
            self.gff_df = pd.concat(
                [self.gff_df, pd.read_csv(gff_path, names=self.column_names, comment="#", sep="\t")],
                ignore_index=True)
        self.gff_df.reset_index(drop=True, inplace=True)
        self.regions = pd.concat([self.regions, self.gff_df[self.gff_df["type"] == "region"]])
        self.gff_df.drop(self.regions.index, inplace=True, axis=0)
        print(f"==> Parsed {len(parsed_paths)} GFFs")

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
