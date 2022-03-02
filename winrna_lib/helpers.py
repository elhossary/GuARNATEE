import sys
import typing
import numpy as np
import pandas as pd
from Bio import SeqUtils


class Helpers:
    def __init__(self):
        pass

    @staticmethod
    def expand_attributes_to_columns(in_df) -> pd.DataFrame:
        df = in_df.copy()
        df.sort_values(["seqid", "start", "end"], inplace=True)
        df.reset_index(inplace=True, drop=True)
        if "attributes" in df.columns:
            df = Helpers.explode_dict_yielding_func_into_columns(df, "attributes", Helpers.parse_attributes)
            """
            if "attributes" in df.columns:
                for i in df.index:
                    attr_str = df.at[i, "attributes"]
                    if attr_str != "":
                        for k, v in Helpers.parse_attributes(attr_str).items():
                            df.at[i, k] = v
            """
            df.drop(["attributes"], inplace=True, axis=1)

        else:
            print("Warning: Attributes column not found!")
        return df

    @staticmethod
    def warp_non_gff_columns(gff_df: pd.DataFrame, keep_columns=False) -> pd.DataFrame:
        gff_df = gff_df.copy()
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
        non_gff_columns = [x for x in gff_df.columns.tolist() if x not in gff_columns]
        if len(non_gff_columns) == 0:
            return gff_df
        for i in gff_df.index:
            extra_attr = [f"{col}={gff_df.at[i, col]}" for col in non_gff_columns if gff_df.at[i, col] not in ["", np.nan, np.NAN, "nan", None]]
            extra_attr = [x for x in extra_attr if not x.endswith("=")]
            #print(extra_attr)
            gff_df.at[i, "extra_attributes"] = ";".join(extra_attr).strip(";")
        if "attributes" in gff_df.columns.tolist():
            gff_df["attributes"] = (
                gff_df["attributes"] + ";" + gff_df["extra_attributes"]
            )
            non_gff_columns.append("extra_attributes")
        else:
            gff_df.rename(columns={"extra_attributes": "attributes"}, inplace=True)
        if not keep_columns:
            gff_df.drop(non_gff_columns, inplace=True, axis=1)
            gff_df = gff_df.reindex(columns=gff_columns)
        gff_df["attributes"] = gff_df["attributes"].str.strip(to_strip=";")
        gff_df["attributes"] = gff_df["attributes"].str.replace(";;", ";")
        return gff_df

    @staticmethod
    def get_gff_df(df, anno_source="", anno_type="", strand="", new_id=False):
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
        essential_columns = ["seqid", "start", "end"]
        df_columns = df.columns.tolist()
        # Checking inputs
        if not set(essential_columns).issubset(df_columns):
            print("Error: Missing essential columns")
            sys.exit(1)
        if "strand" not in df_columns and strand not in ["+", "-"]:
            print("Error: Missing strand information")
            sys.exit(1)
        # Adding missing columns
        if "source" not in df_columns:
            df["source"] = anno_source
        if "type" not in df_columns:
            df["type"] = anno_type
        if "score" not in df_columns:
            df["score"] = "."
        if "phase" not in df_columns:
            df["phase"] = "."
        if "strand" not in df_columns:
            df["strand"] = strand
        if "attributes" not in df_columns:
            df["attributes"] = ""
        df = df.reindex(columns=gff_columns)

        for i in df.index:
            strand_letter = "F" if df.at[i, "strand"] == "+" else "R"
            if new_id:
                old_attr_str = df.at[i, "attributes"]
                attr = (
                    Helpers.parse_attributes(old_attr_str) if old_attr_str != "" else {}
                )
                if "id" in attr.keys():
                    del attr["id"]
                if "name" in attr.keys():
                    del attr["name"]
                df.at[i, "attributes"] = (
                    f'ID={df.at[i, "seqid"]}{strand_letter}_{anno_type}_{i}'
                    f';name={df.at[i, "seqid"]}{strand_letter}_{anno_type}_{i}'
                )
                attr_addition = ";".join([f"{k}={v}" for k, v in attr.items()])
                if attr_addition != "":
                    df.at[i, "attributes"] += ";" + attr_addition
        df["attributes"] = df["attributes"].str.strip(to_strip=";")
        df = Helpers.warp_non_gff_columns(df)
        df.sort_values(["seqid", "start", "end"], inplace=True)
        return df

    @staticmethod
    def get_naming_attributes(attr_str: str, prefix="") -> str:
        accepted_attributes = ["name", "gene", "locus_tag", "old_locus_tag"]
        attr_dict = Helpers.parse_attributes(attr_str)
        drop_ids = []
        for k in attr_dict.keys():
            if not any(aa in k for aa in accepted_attributes):
                drop_ids.append(k)
        for k in drop_ids:
            del attr_dict[k]
        return ";".join([f"{prefix}{k}={v}" for k, v in attr_dict.items()])

    @staticmethod
    def parse_attributes(attr_str: str):
        if attr_str == "":
            return {}
        attr_pairs = attr_str.split(";")
        attr_dict = {}
        for attr_pair in attr_pairs:
            attr_pair_lst = attr_pair.split("=")
            if len(attr_pair_lst) != 2:
                print(f"Warning: Skipping ambiguous key/value pair in GFF at: {attr_str}")
                continue
            k, v = attr_pair_lst[0], attr_pair_lst[1]
            if v in attr_dict.values():
                continue
            if k.lower() in attr_dict.keys():
                if attr_dict[k.lower()] == v:
                    continue
                attr_dict[k.lower()] += "," + v
            else:
                attr_dict[k.lower()] = v
        return attr_dict

    @staticmethod
    def flatten_attr_dict(in_dict):
        return ";".join([f"{k}={v}" for k, v in in_dict.items()])

    @staticmethod
    def get_gc_content(seq_str):
        return SeqUtils.GC(seq_str)

    @staticmethod
    def explode_dict_yielding_func_into_columns(df: pd.DataFrame, df_col: str, func: typing.Callable, column_prefix="") -> pd.DataFrame:
        df.reset_index(inplace=True, drop=True)
        df["TMP_COLUMN"] = df[df_col].map(lambda x: func(x))
        tmp_df = df["TMP_COLUMN"].apply(pd.Series)
        tmp_df.fillna("", inplace=True)
        if column_prefix != "":
            tmp_df = tmp_df.add_prefix(column_prefix)
        df = pd.merge(left=df, right=tmp_df, right_index=True, left_index=True, how='left')
        df.drop(columns=["TMP_COLUMN"], inplace=True)
        return df
