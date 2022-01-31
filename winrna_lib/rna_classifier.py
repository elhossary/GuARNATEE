import pandas as pd
import pybedtools as pybed
from itertools import product
from io import StringIO
from tqdm import tqdm

class RNAClassifier:

    def __init__(self, gff_df, anno_tbl_df):
        self.gff_df = gff_df[gff_df["type"] != "region"]
        self.anno_tbl_df = anno_tbl_df
        self.gff_columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
        self.classes = pd.DataFrame(columns=self.gff_columns)
        self.classify()

    def classify(self):
        self._classify()
        self._drop_redundancies()
        self.classes["source"] = "WinRNA"
        self.classes.sort_values(["seqid", "start", "end"], inplace=True)
        #print(self.classes.to_string())

    def _drop_redundancies(self):
        df = self.expand_attributes_to_columns(self.classes)
        combs = list(product(df["seqid"].unique().tolist(), ["+", "-"], ["start", "end"]))
        drop_ids = []
        for seqid, strand, select_column in tqdm(combs, desc="===> Cleaning redundancies"):
            df_slice = df[(df["seqid"] == seqid) & (df["strand"] == strand)]
            select_keys = df_slice[select_column].unique().tolist()
            for select_key in select_keys:
                tmp_df = df_slice[df_slice[select_column] == select_key].copy()
                if tmp_df.shape[0] == 1:
                    continue
                if tmp_df["type"].unique().size > 0:
                    # Split
                    print(tmp_df.to_string())
                tmp_df.sort_values(by=["SS_step_factor", "TS_step_factor", "length", "SS_height", "TS_height"],
                                   ascending=[False, False, True, False, False],
                                   inplace=True)
                drop_ids.extend(tmp_df.index.tolist()[1:])

        self.classes.drop(drop_ids, axis=0, inplace=True)

    @staticmethod
    def expand_attributes_to_columns(df) -> pd.DataFrame:
        df.sort_values(["seqid", "start", "end"], inplace=True)
        df.reset_index(inplace=True, drop=True)
        for i in df.index:
            attributes = df.at[i, "attributes"].split(";")
            for attr in attributes:
                k, v = attr.split("=")
                df.at[i, k] = v
        df.drop(["attributes"], inplace=True, axis=1)
        return df

    def _classify(self) -> None:
        # full_ref_df = self.gff_df[~self.gff_df["type"].isin(["CDS", "exon", "protein_binding_site"])].copy()
        # [self.gff_df["type"] == "gene"].copy()
        full_ref_df = self.gff_df
        anno_tbl_bed = pybed.BedTool.from_dataframe(self.anno_tbl_df)
        # drop unwanted overlaps like tRNA or rRNA
        drop_ref_df = full_ref_df[(full_ref_df["type"] == "tRNA") |
                                  (full_ref_df["attributes"].str.contains("gene_biotype=tRNA;")) |
                                  (full_ref_df["type"] == "rRNA") |
                                  (full_ref_df["attributes"].str.contains("gene_biotype=rRNA;"))].copy()
        drop_ref_bed = pybed.BedTool.from_dataframe(drop_ref_df)
        anno_tbl_bed = anno_tbl_bed.intersect(drop_ref_bed, v=True, f=0.10, s=True)
        #
        ncrna_ref_df = full_ref_df[full_ref_df["type"] == "ncRNA"].copy()
        cds_genes_ref_df = full_ref_df[(full_ref_df["type"] == "gene") &
                                       (full_ref_df["attributes"].str.contains("gene_biotype=protein_coding;"))].copy()

        # get novel intergenic (non-overlaps of any other annotations)
        full_ref_bed = pybed.BedTool.from_dataframe(full_ref_df)
        intergenic_df = \
            anno_tbl_bed.intersect(full_ref_bed, v=True, f=0.99, s=True).to_dataframe(names=self.gff_columns)
        anno_tbl_bed = pybed.BedTool.from_dataframe(self.subtract_dfs(self.anno_tbl_df, intergenic_df))
        intergenic_df["type"] = "intergenic_srna_candidate"
        intergenic_df["attributes"] += ";annotation_type=intergenic_sRNA_candidate"
        self.classes = pd.concat([self.classes, intergenic_df], ignore_index=True)

        # get known ncRNAs
        if ncrna_ref_df.shape[0] > 0:
            ncrna_ref_bed = pybed.BedTool.from_dataframe(ncrna_ref_df.loc[:, self.gff_columns[:-1]])
            ncrna_bed = anno_tbl_bed.intersect(ncrna_ref_bed, wa=True, f=0.75, u=True, r=True, s=True)
            ncrna_df = self.pybed_to_df_func(ncrna_bed, self.gff_columns)
            ncrna_df = self.map_annotations(ncrna_ref_df, ncrna_df, "intergenic")
            self.classes = pd.concat([self.classes, ncrna_df], ignore_index=True)

        # get unknown ORF_int
        if cds_genes_ref_df.shape[0] > 0:
            cds_genes_ref_bed = pybed.BedTool.from_dataframe(cds_genes_ref_df.loc[:, self.gff_columns[:-1]])
            orf_int_bed = anno_tbl_bed.intersect(cds_genes_ref_bed, wa=True, f=0.50, u=True, s=True)
            orf_int_df = self.pybed_to_df_func(orf_int_bed, self.gff_columns)
            orf_int_df = self.map_annotations(cds_genes_ref_df, orf_int_df, "orf_int")
            self.classes = pd.concat([self.classes, orf_int_df], ignore_index=True)


    @staticmethod
    def subtract_dfs(full_df: pd.DataFrame, partial_df: pd.DataFrame):
        return pd.merge(left=partial_df, right=full_df, indicator=True, how='outer')\
            .query('_merge=="right_only"').drop('_merge', axis=1)

    def map_annotations(self, ref_df: pd.DataFrame, df: pd.DataFrame, mode: str):
        ref_df["source"] = "ref"
        df["source"] = "query"
        all_df = pd.concat([ref_df, df], ignore_index=True)
        all_bed = pybed.BedTool.from_dataframe(all_df)
        all_bed = all_bed.cluster(d=0, s=True)
        all_df = self.pybed_to_df_func(all_bed, self.gff_columns + ["cluster"])
        ref_df = all_df[all_df["source"] == "ref"]
        query_df = all_df[all_df["source"] == "query"]
        query_df = pd.merge(left=query_df, right=ref_df.loc[:, ["cluster", "type", "start", "end", "attributes"]],
                            on="cluster", how='left', suffixes=("", "_ref"))
        if mode == "intergenic":
            query_df = self.mark_known_intergenic(query_df)
        elif mode == "orf_int":
            query_df = self.mark_known_orf_int(query_df)
        else:
            print("internal error")
            exit(1)
        query_df.drop(columns=["cluster", "type_ref", "start_ref", "end_ref", "attributes_ref"], inplace=True)
        return query_df

    def mark_known_orf_int(self, df: pd.DataFrame):
        df.reset_index(inplace=True, drop=True)
        drop_ids = []
        for i in df.index:
            check_range = range(df.at[i, "start_ref"], df.at[i, "end_ref"] + 1, 1)
            if df.at[i, "start"] in check_range and df.at[i, "end"] in check_range:
                ref_attr = self.parse_attributes(df.at[i, "attributes_ref"])
                locus_tag = ref_attr["locus_tag"] if "locus_tag" in ref_attr.keys() else "NA"
                gene = ref_attr["name"] if "name" in ref_attr.keys() else "NA"
                #gene = ref_attr["gene"] if "gene" in ref_attr.keys() and gene == "NA" else "NA"
                df.at[i, "type"] = "ORF_int_candidate"
                df.at[i, "attributes"] += \
                    f";annotation_type=ORF_internal_sRNA;host_gene_name={gene};sRNA_locus_tag={locus_tag}"
            else:
                drop_ids.append(i)
        df.drop(index=drop_ids, inplace=True)
        return df

    def mark_known_intergenic(self, df: pd.DataFrame):
        for i in df.index:
            ref_attr = self.parse_attributes(df.at[i, "attributes_ref"])
            locus_tag = ref_attr["locus_tag"] if "locus_tag" in ref_attr.keys() else "NA"
            gene = ref_attr["gene"] if "gene" in ref_attr.keys() else "NA"
            df.at[i, "type"] = "ncRNA"
            df.at[i, "attributes"] += \
                f";annotation_type=known_intergenic_sRNA;sRNA_name={gene};sRNA_locus_tag={locus_tag}"
        return df

    @staticmethod
    def pybed_to_df_func(pybed_obj: pybed, names: list) -> pd.DataFrame:
        gff_columns = {"seqname": "seqid", "feature": "type", "frame": "phase"}
        ret_df = pybed_obj.to_dataframe(names=names)
        ret_df.rename(columns=gff_columns, inplace=True)
        return ret_df

    @staticmethod
    def parse_attributes(attr_str):
        return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}