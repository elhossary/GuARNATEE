import pandas as pd
import numpy as np


class Helpers:

    def __init__(self):
        pass

    @staticmethod
    def expand_attributes_to_columns(df) -> pd.DataFrame:
        df.sort_values(["seqid", "start", "end"], inplace=True)
        df.reset_index(inplace=True, drop=True)
        if "attributes" in df.columns:
            for i in df.index:
                attr_str = df.at[i, "attributes"]
                if attr_str != "":
                    for k, v in Helpers.parse_attributes(attr_str).items():
                        df.at[i, k] = v
            df.drop(["attributes"], inplace=True, axis=1)
        else:
            print("Warning: Attributes column not found!")
        return df

    @staticmethod
    def parse_attributes(attr_str):
        attr_str = attr_str.strip(";")
        return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}

    @staticmethod
    def warp_non_gff_columns(gff_df: pd.DataFrame) -> pd.DataFrame:
        gff_columns = \
            ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
        non_gff_columns = [x for x in gff_df.columns.tolist() if x not in gff_columns]
        if len(non_gff_columns) == 0:
            return gff_df
        for i in gff_df.index:
            extra_attr = [f'{col}={gff_df.at[i, col]}' for col in non_gff_columns]
            gff_df.at[i, "extra_attributes"] = \
                ";".join([x for x in extra_attr if "=nan" not in x]).strip(";")
        if "attributes" in gff_df.columns.tolist():
            gff_df["attributes"] = gff_df["attributes"] + ";" + gff_df["extra_attributes"]
            non_gff_columns.append("extra_attributes")
        else:
            gff_df.rename(columns={"extra_attributes": "attributes"}, inplace=True)

        gff_df.drop(non_gff_columns, inplace=True, axis=1)
        gff_df = gff_df.reindex(columns=gff_columns)
        gff_df["attributes"] = gff_df["attributes"].str.strip(to_strip=";")
        return gff_df

    @staticmethod
    def get_gff_df(df, anno_source="", anno_type="", strand="", new_id=False):
        gff_columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
        essential_columns = ["seqid", "start", "end"]
        df_columns = df.columns.tolist()
        # Checking inputs
        if not set(essential_columns).issubset(df_columns):
            print("Error: Missing essential columns")
            exit(1)
        if "strand" not in df_columns and strand not in ["+", "-"]:
            print("Error: Missing strand information")
            exit(1)
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
                attr = Helpers.parse_attributes(old_attr_str) if old_attr_str != "" else {}
                if "id" in attr.keys():
                    del attr["id"]
                if "name" in attr.keys():
                    del attr["name"]
                df.at[i, "attributes"] = f'ID={df.at[i, "seqid"]}{strand_letter}_{anno_type}_{i}' \
                                         f';name={df.at[i, "seqid"]}{strand_letter}_{anno_type}_{i}'
                attr_addition = ";".join([f'{k}={v}' for k, v in attr.items()])
                if attr_addition != "":
                    df.at[i, "attributes"] += ";" + attr_addition
        df["attributes"] = df["attributes"].str.strip(to_strip=";")
        df = Helpers.warp_non_gff_columns(df)
        df.sort_values(["seqid", "start", "end"], inplace=True)
        return df
