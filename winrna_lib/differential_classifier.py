import pandas as pd
import pybedtools as pybed
import numpy as np
from winrna_lib.helpers import Helpers


class DifferentialClassifier:
    """
    This class is to compare two sets of prediction typically TEX+ and TEX-
    """

    def __init__(self, in_dict: dict):
        # in_dict takes a pair of dataframes in as values of a dict where the keys are the lib names
        self.in_dict = in_dict
        if len(in_dict.keys()) != 2:
            print("Error: non-uniformed sets passed")

    def score_similarity(self):
        dfs = list(self.in_dict.items())
        df1_lib, df2_lib = dfs[0][0], dfs[1][0]
        df1, df2 = dfs[0][1], dfs[1][1]
        bed1_score_df = self._score_similarity(df1, df2, df1_lib, df2_lib)
        bed2_score_df = self._score_similarity(df2, df1, df1_lib, df1_lib)
        return bed1_score_df, bed2_score_df

    def _score_similarity(
        self, df1: pd.DataFrame, df2: pd.DataFrame, df1_lib_name: str, df2_lib_name: str
    ) -> pd.DataFrame:
        """ """
        gff_columns = [
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
        scores_columns = [
            "overlaps_count",
            "overlap_size",
            "length",
            "overlap_percentage",
        ]
        for i in df1.index:
            df1.at[i, "attributes"] += f";lib_type={df1_lib_name}"
        bed1 = pybed.BedTool.from_dataframe(df1).sort()
        bed2 = pybed.BedTool.from_dataframe(df2).sort()
        bed_score_df = bed1.coverage(bed2, s=True).to_dataframe(
            names=gff_columns + scores_columns
        )
        bed_score_df.drop(["overlap_size", "length"], inplace=True, axis=1)
        bed_score_df["overlap_percentage"] = bed_score_df["overlap_percentage"] * 100
        bed_score_df["overlap_percentage"] = bed_score_df["overlap_percentage"].round(2)
        bed_score_df["overlap_percentage"].replace(
            to_replace=0.0, value=np.nan, inplace=True
        )
        bed_score_df["overlaps_count"] = bed_score_df["overlaps_count"].astype(int)
        # bed_score_df["overlaps_count"].replace(to_replace=0.0, value=np.nan, inplace=True)
        bed_score_df.rename(
            columns={
                "overlap_percentage": f"{df2_lib_name}_overlap_percentage",
                "overlaps_count": f"{df2_lib_name}_overlaps_count",
            },
            inplace=True,
        )
        bed_score_df = Helpers.warp_non_gff_columns(bed_score_df)
        return bed_score_df
