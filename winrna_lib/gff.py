import glob
import os.path

import numpy as np
import pandas as pd
from winrna_lib.helpers import Helpers


class GFF:
    def __init__(self, gff_paths: list):
        self.gff_paths = gff_paths
        self.column_names = [
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
        self.gff_df = pd.DataFrame(columns=self.column_names)
        self.regions = pd.DataFrame()
        self.parse()
        self.seqid_groups = {}

    def parse(self):
        print("=> Parsing input GFF file")
        parsed_paths = []
        for item in self.gff_paths:
            for sub_item in glob.glob(item):
                parsed_paths.append(os.path.abspath(sub_item))
        for gff_path in parsed_paths:
            gff_parsed = pd.read_csv(
                gff_path, names=self.column_names, comment="#", sep="\t"
            )
            # self.seqid_groups[gff_path] = gff_parsed["seqid"].unique().tolist()
            self.gff_df = pd.concat([self.gff_df, gff_parsed], ignore_index=True)
        self.gff_df.reset_index(drop=True, inplace=True)
        self.regions = pd.concat(
            [self.regions, self.gff_df[self.gff_df["type"] == "region"]]
        )
        self.gff_df.drop(self.regions.index, inplace=True, axis=0)
        print(f"==> Parsed {len(parsed_paths)} GFF files")

