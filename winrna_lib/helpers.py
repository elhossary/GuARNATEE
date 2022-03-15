import sys
import typing
import numpy as np
import pandas as pd
from Bio import SeqUtils
from multiprocessing import Pool


class Helpers:
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
    def warp_non_gff_columns(gff_df: pd.DataFrame, exclude_columns=None, no_join=False, keep_columns=False) -> pd.DataFrame:
        gff_df = gff_df.copy()
        gff_df.reset_index(drop=True, inplace=True)
        if exclude_columns is None:
            exclude_columns = []
        gff_columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
        non_gff_columns = [x for x in gff_df.columns.tolist() if x not in gff_columns + exclude_columns]
        if len(non_gff_columns) == 0:
            return gff_df
        #prohibited_values = ["", np.nan, np.NAN, "nan", None]
        #attr_dict_list = gff_df.loc[:, non_gff_columns].to_dict('records')
        for i in gff_df.index:
            gff_df.at[i, "extra_attributes"] = Helpers.attributes_dict_to_str({col: gff_df.at[i, col] for col in non_gff_columns})
        #gff_df["extra_attributes"] = \
        #    pd.Series(list(map(Helpers.attributes_dict_to_str, attr_dict_list)))
        #print(gff_df[gff_df["start"] == 883439].to_string())
        if not no_join and "extra_attributes" in gff_df.columns:
            gff_df['attributes'] = gff_df['attributes'] + ";" + gff_df['extra_attributes']
            non_gff_columns.append("extra_attributes")
        if not keep_columns:
            gff_df.drop(non_gff_columns, inplace=True, axis=1)
            gff_df = gff_df.reindex(columns=gff_columns)
        gff_df["attributes"] = gff_df["attributes"].str.strip(to_strip=";")
        return gff_df

    @staticmethod
    def attributes_dict_to_str(attr_dict):
        prohibited_values = ["", None, "nan", np.nan]
        return ";".join([f"{k}={v}" for k, v in attr_dict.items() if v not in prohibited_values])

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
            sys.sys.exit(1)
        if "strand" not in df_columns and strand not in ["+", "-"]:
            print("Error: Missing strand information")
            sys.sys.exit(1)
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
        seq_lst = df[df_col].tolist()
        df["TMP_COLUMN"] = pd.Series(list(Pool().map(func, seq_lst)))
        return Helpers.explode_column_of_dicts(df, "TMP_COLUMN", column_prefix)

    @staticmethod
    def explode_column_of_dicts(df: pd.DataFrame, df_col: str, column_prefix="") -> pd.DataFrame:
        df.reset_index(inplace=True, drop=True)
        tmp_df = df[df_col].apply(pd.Series)
        tmp_df.fillna("", inplace=True)
        if column_prefix != "":
            tmp_df = tmp_df.add_prefix(column_prefix)
        df = pd.merge(left=df, right=tmp_df, right_index=True, left_index=True, how='left')
        df.drop(columns=["TMP_COLUMN"], inplace=True)
        return df

    @staticmethod
    def parse_attributes_into_dict(gff_df: pd.DataFrame, attr_col="attributes") -> pd.DataFrame:
        gff_df[f"{attr_col}_dict"] = gff_df[attr_col].map(Helpers.parse_attributes_str)
        return gff_df

    @staticmethod
    def add_type_as_prefix_to_attributes_keys(gff_df: pd.DataFrame, attr_col="attributes", type_col="type") -> pd.DataFrame:
        gff_df.reset_index(inplace=True, drop=True)
        gff_df = Helpers.parse_attributes_into_dict(gff_df, attr_col)
        for i in gff_df.index:
            attr_type = gff_df.at[i, type_col]
            old_dict = gff_df.at[i, f"{attr_col}_dict"]
            new_dict = {(f"{attr_type}_{k}" if attr_type.lower() not in k.lower() else k): v for k, v in old_dict.items()}
            gff_df.at[i, attr_col] = Helpers.attributes_dict_to_str(new_dict)
        gff_df.drop(columns=[f"{attr_col}_dict"], inplace=True)
        return gff_df

    @staticmethod
    def merge_same_intervals(gff_df: pd.DataFrame):
        column_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
        gff_df = Helpers.add_type_as_prefix_to_attributes_keys(gff_df)
        essential_columns = ["seqid", "start", "end", "strand"]
        gff_df = gff_df.groupby(essential_columns, as_index=False).agg({"source": list,
                                                                        "type": list,
                                                                        "phase": list,
                                                                        "score": list,
                                                                        "attributes": list})
        for i in gff_df.index:
            for column in ["source", "type", "phase", "score"]:
                gff_df.at[i, column] = "|".join([str(x) for x in set(gff_df.at[i, column])])
            gff_df.at[i, "attributes"] = ";".join(set(gff_df.at[i, "attributes"]))
        gff_df.reset_index(inplace=True, drop=True)
        gff_df = gff_df.reindex(columns=column_names)
        return gff_df

    @staticmethod
    def parse_attributes_str(attr_str):
        if attr_str in ["", None]:
            return {}
        attr_pairs = [attr.split("=") for attr in attr_str.split(";") if attr != ""]
        attr_dict = {}
        for attr_pair in attr_pairs:
            if len(attr_pair) != 2:
                print(f"Warning: Skipping ambiguous key/value pair in GFF at: {attr_pairs}")
                continue
            if attr_pair[0] in attr_dict.keys():
                if attr_pair[1] == attr_dict[attr_pair[0]]:
                    continue
                attr_dict[attr_pair[0]] += f"|{attr_pair[1]}"
                continue
            attr_dict[attr_pair[0]] = attr_pair[1]
        return attr_dict

    @staticmethod
    def filter_attributes(gff_df: pd.DataFrame, filters: list, attr_col="attributes") -> pd.DataFrame:
        gff_df = Helpers.parse_attributes_into_dict(gff_df, attr_col)
        for i in gff_df.index:
            gff_df.at[i, f"{attr_col}_dict"] = \
                {k: v for k, v in gff_df.at[i, f"{attr_col}_dict"].items() if
                 any(f.lower() in k.lower() for f in filters)}
        gff_df[attr_col] = gff_df[f"{attr_col}_dict"].apply(Helpers.attributes_dict_to_str)
        gff_df.drop(columns=[f"{attr_col}_dict"], inplace=True)
        return gff_df

    @staticmethod
    def rewrap_attributes_column(gff_df, attr_col="attributes") -> pd.DataFrame:
        gff_df = Helpers.parse_attributes_into_dict(gff_df, attr_col)
        for i in gff_df.index:
            gff_df.at[i, attr_col] = Helpers.attributes_dict_to_str(gff_df.at[i, f"{attr_col}_dict"])
        gff_df.drop(columns=[f"{attr_col}_dict"], inplace=True)
        return gff_df
