import pandas as pd
import pybedtools as pb
from tqdm import tqdm
from winrna_lib.helpers import Helpers


class CandidatesMerger:

    def __init__(self, in_df: pd.DataFrame, similarity=0.85):
        self.column_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
        self.in_df = self.cluster_similar_annotations(in_df, similarity)
        self.in_df = self.merge_by_group_by(self.in_df)

    def merge(self):
        return self.in_df

    def merge_by_group_by(self, df: pd.DataFrame):
        df = Helpers.expand_attributes_to_columns(df)
        #df["condition_full_name"] = df["condition"] + "_" + df["lib_type"]
        df["merge_count"] = "none"
        ret_df = pd.DataFrame()
        merge_cols = ["gene_name", "gene_locus_tag", "cds_protein_id",
                      "annotation_class", "sub_class", "detection_status", "ncrna_name", "cluster_id"]
        num_cols = ["sum_all_rank", "sum_step_factors_mfe_rank", "mfe_rank", "step_factors_rank",
                    "plateau_heights_rank", "ss_step_factor", "ts_step_factor", "mfe"]
        for col in num_cols:
            df[col] = pd.to_numeric(df[col], errors='coerce', downcast='float')
        for group in tqdm(df["group_id"].unique(), desc="Merging clusters"):
            tmp_df = df[df["group_id"] == group].copy()
            ret_df = pd.concat([ret_df, tmp_df.groupby(self.column_names[:-1] + merge_cols, as_index=False)
                               .agg({'start': min,
                                     'end': max,
                                     "sum_all_rank": max,
                                     "sum_step_factors_mfe_rank": max,
                                     "merge_count": 'count',
                                     "mfe_rank": max,
                                     "mfe": min,
                                     "step_factors_rank": max,
                                     "ss_step_factor": max,
                                     "ts_step_factor": max,
                                     "plateau_heights_rank": max})], ignore_index=True)
        ret_df["id"] = [f"cand_{x}" for x in range(1, ret_df.shape[0] + 1, 1)]
        ret_df["name"] = [f"cand_{x}" for x in range(1, ret_df.shape[0] + 1, 1)]
        ret_df = Helpers.warp_non_gff_columns(ret_df)
        return ret_df

    def cluster_similar_annotations(self, in_df, similarity, passes=2, pass_num=0):
        if "group_id" in in_df.columns:
            in_df.drop(columns=["group_id"], inplace=True)
        ret_df = pd.DataFrame()
        counter = 0
        skip_indexes = []
        in_pb = pb.BedTool().from_dataframe(in_df).sort()
        pass_num += 1
        for i in tqdm(in_df.index, desc=f'Clustering similar intervals, pass: {pass_num}'):
            if i in skip_indexes:
                continue
            irow = in_df.loc[[i]]

            similar_intervals_pb = in_pb.intersect(pb.BedTool().from_dataframe(irow), s=True, f=similarity, r=True, wa=True)

            similar_intervals_df = similar_intervals_pb.to_dataframe(names=self.column_names)
            if similar_intervals_df.empty:
                continue
            counter += 1
            similar_intervals_df["group_id"] = counter
            ret_df = pd.concat([ret_df, similar_intervals_df], ignore_index=True)
            similar_intervals_df.drop(columns=["group_id"], inplace=True)
            skip_indexes.extend(
                pd.merge(left=in_df, right=similar_intervals_df, on=self.column_names, how='left', indicator=True)
                    .query("_merge=='both'").index.tolist())

        if pass_num < passes:
            ret_df = self.cluster_similar_annotations(ret_df, similarity, passes, pass_num)
        return ret_df