import sys
import os
import io
from itertools import product
import tempfile
from Bio import SeqUtils
import numpy as np
import pandas as pd
import pybedtools as pybed
import RNA
from tqdm import tqdm
from more_itertools import consecutive_groups
from winrna_lib.helpers import Helpers
from winrna_lib.fasta import Fasta


class RNAClassifier:
    def __init__(self, gff_df, anno_tbl_df, fasta: Fasta):
        self.gff_df = gff_df[gff_df["type"] != "region"].copy()
        self.anno_tbl_df = anno_tbl_df
        self.fasta = fasta
        self.gff_columns = [
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
        self.seqids = set.intersection(set(self.gff_df["seqid"].unique().tolist()),
                                       set(self.fasta.fwd_seqs.keys()),
                                       set(self.anno_tbl_df["seqid"].unique().tolist()))
        self.anno_tbl_df = self.anno_tbl_df[self.anno_tbl_df["seqid"].isin(self.seqids)]
        self.classes = pd.DataFrame
        self.dispatch()  # run the classifier at init

    def dispatch(self):
        self.prefilter_candidates()
        self.classify()
        self.get_intergenic_flanks()
        self.get_gff_sequences_features(is_rna=True)
        #print(self.classes.head(20).to_string())
        self._drop_redundancies()
        self.classes.sort_values(["seqid", "start", "end"], inplace=True)

    def _drop_redundancies(self):
        df = Helpers.expand_attributes_to_columns(self.classes)
        combs = list(
            product(df["seqid"].unique().tolist(), ["+", "-"], ["start", "end"])
        )
        drop_ids = []
        for seqid, strand, select_column in tqdm(
            combs, desc="===> Cleaning redundancies"
        ):
            df_slice = df[(df["seqid"] == seqid) & (df["strand"] == strand)]
            select_keys = df_slice[select_column].unique().tolist()
            for select_key in select_keys:
                tmp_df = df_slice[df_slice[select_column] == select_key].copy()
                if tmp_df.shape[0] == 1:
                    continue
                # TODO
                """
                if tmp_df["type"].unique().size > 1:
                    # Split
                    print(tmp_df.to_string())
                    continue
                """
                tmp_df.sort_values(
                    by=[
                        "ss_step_factor",
                        "ts_step_factor",
                        "length",
                        "ss_height",
                        "ts_height",
                    ],
                    ascending=[False, False, True, False, False],
                    inplace=True,
                )
                drop_ids.extend(tmp_df.index.tolist()[1:])
        self.classes.drop(drop_ids, axis=0, inplace=True)
        # self._drop_redundacies_in_clusters()

    def _drop_redundacies_in_clusters(self):
        # TODO
        clusteres_df = (
            pybed.BedTool.from_dataframe(self.classes)
            .sort()
            .cluster(s=True, d=0)
            .to_dataframe(names=self.gff_columns + ["cluster"])
        )
        # print(clusteres_df.to_string())

    def classify(self):
        print("=> Classifying candidates")
        classes_df = self.anno_tbl_df.copy()
        classes_df["type"] = "sRNA_candidate"
        classes_df["annotation_class"] = ""
        classes_df["status"] = ""

        ref_df = self.re_classify_ref_gff()
        ref_pb = pybed.BedTool.from_dataframe(ref_df).sort()
        classes_pb = pybed.BedTool.from_dataframe(
            Helpers.warp_non_gff_columns(classes_df)
        ).sort()
        ref_columns = (
            self.gff_columns
            + [f"ref_{col}" for col in self.gff_columns]
            + ["overlap_size"]
        )

        classes_df = (
            classes_pb.intersect(ref_pb, wao=True, f=0.5)
            .to_dataframe(names=ref_columns)
            .drop(columns=["ref_seqid", "ref_source", "ref_score", "ref_phase"])
        )

        ref_columns = [c for c in classes_df.columns if "ref_" in c]
        for rc in ref_columns:
            classes_df.loc[classes_df[rc] == ".", rc] = ""
            #classes_df.loc[classes_df[rc] == -1, rc] = 0

        classes_df = classes_df.groupby(self.gff_columns + ["ref_strand"], as_index=False).agg({
            "overlap_size": max, "ref_attributes": ";".join, "ref_type": "_".join, "ref_end": max, "ref_start": min})
        classes_df.replace("ncRNA_ncRNA", "ncRNA", inplace=True)
        #classes_df = classes_df.loc[classes_df.groupby(self.gff_columns).overlap_size.idxmax()]
        for i in classes_df.index:
            classes_df.at[i, "ref_attributes"] = Helpers.get_naming_attributes(classes_df.at[i, "ref_attributes"])
        classes_df = classes_df.reindex(columns=self.gff_columns + [c for c in classes_df.columns if c not in self.gff_columns])
        classes_df.loc[
            classes_df["overlap_size"] == 0, ["annotation_class", "status"]
        ] = ("intergenic", "novel")
        classes_df.loc[
            (classes_df["ref_type"] == "CDS")
            & (classes_df["strand"] == classes_df["ref_strand"]),
            ["annotation_class", "status"],
        ] = ("ORF_int", "novel")

        classes_df.loc[
            (classes_df["ref_type"] == "CDS")
            & (classes_df["strand"] != classes_df["ref_strand"]),
            ["annotation_class", "status"]
        ] = ("antisense_to_CDS_region", "novel")

        classes_df.loc[
            (classes_df["ref_type"].str.contains("ncRNA"))
            & (classes_df["strand"] == classes_df["ref_strand"]),
            ["annotation_class", "status"]
        ] = ("ncRNA", "known")
        classes_df.loc[
            (classes_df["ref_type"].str.contains("ncRNA"))
            & (classes_df["strand"] != classes_df["ref_strand"]),
            ["annotation_class", "status"]
        ] = ("antisense_to_ncRNA_region", "novel")

        classes_df.loc[
            (classes_df["ref_type"].str.contains("ORF_int"))
            & (classes_df["strand"] == classes_df["ref_strand"]),
            ["annotation_class", "status", "ref_type"]
        ] = ("ORF_int", "known", "CDS")

        classes_df = classes_df.apply(self.add_overlap_info, axis=1)
        classes_df.loc[
            (classes_df["ref_type"] == "CDS")
            & (classes_df["annotation_class"] == "ORF_int")
            & (classes_df["overlap_fragment_ratio"] >= 75),
            ["annotation_class", "status"]
        ] = ("sORF", "known")

        classes_df["annotation_class"].fillna("bbNA", inplace=True)
        classes_df["status"].fillna("NA", inplace=True)
        classes_df["attributes"] = classes_df["attributes"].str.cat(classes_df["gene_info"].fillna(""), sep=";")
        classes_df.drop(
            columns=["gene_info", "ref_attributes", "overlap_size", "ref_type", "ref_end", "ref_start", "ref_strand"],
            inplace=True
        )
        classes_df.fillna("", inplace=True)
        classes_df.reset_index(inplace=True, drop=True)
        classes_df = Helpers.warp_non_gff_columns(classes_df)
        classes_df.sort_values(["seqid", "start", "end"], inplace=True)
        self.classes = classes_df

    def re_classify_ref_gff(self):
        cds_df = self.gff_df[self.gff_df["type"] == "CDS"].copy()
        ncrna_df = self.gff_df[self.gff_df["type"] == "ncRNA"].copy()
        cds_pb = pybed.BedTool.from_dataframe(cds_df).sort()
        ncrna_pb = pybed.BedTool.from_dataframe(ncrna_df).sort()
        orf_int_df = ncrna_pb.intersect(cds_pb, wa=True, s=True, f=0.50).to_dataframe(
            names=self.gff_columns
        )
        orf_int_df["type"] = "ORF_int"
        diff_cols = [
            "seqid",
            "source",
            "start",
            "end",
            "score",
            "strand",
            "phase",
            "attributes",
        ]
        ncrna_df = pd.merge(
            left=ncrna_df,
            right=orf_int_df,
            on=diff_cols,
            how="left",
            suffixes=("", "_new"),
        )
        ncrna_df["type"] = ncrna_df["type_new"]
        ncrna_df["type"].fillna("ncRNA", inplace=True)
        ncrna_df.drop(columns=["type_new"], inplace=True)
        return pd.concat([cds_df, ncrna_df], ignore_index=True)

    def add_overlap_info(self, row: pd.Series) -> pd.Series:
        if row["overlap_size"] == 0:
            return row

        overlap_info = self._get_overlap_position(
            (row["start"], row["end"]),
            (row["ref_start"], row["ref_end"]),
        )
        row["overlap_fragment_ratio"] = overlap_info["intersect_size_perc"]
        if row["strand"] == "-":
            overlap_info["diff_before_size_perc"], overlap_info["diff_after_size_perc"] = \
                overlap_info["diff_after_size_perc"], overlap_info["diff_before_size_perc"]
        row["upstream_fragment_ratio"] = overlap_info["diff_before_size_perc"]
        row["downstream_fragment_ratio"] = overlap_info["diff_after_size_perc"]
        row["gene_info"] = Helpers.get_naming_attributes(
            row["ref_attributes"], prefix=f"{row['ref_type']}_"
        )
        return row

    @staticmethod
    def _get_overlap_position(partial: tuple, full: tuple) -> dict:
        partial_set = set(range(partial[0], partial[1] + 1, 1))
        full_set = set(range(full[0], full[1] + 1, 1))
        full_size = len(full_set)
        diff = sorted(full_set.difference(partial_set))
        ret_dict = {
            "intersect_size": 0,
            "diff_before_size": full_size,
            "diff_after_size": 0,
            "intersect_size_perc": 0,
            "diff_before_size_perc": 100,
            "diff_after_size_perc": 0
        }
        if not diff:
            return ret_dict
        intersect = sorted(partial_set.intersection(full_set))
        cga = [list(cg) for cg in consecutive_groups(diff)]  # MUST be a list of 2 lists
        cga_len = len(cga)
        diff_before = []
        diff_after = []
        if cga_len == 0:
            return ret_dict
        elif cga_len == 1:  # Special case
            if max(cga[0]) < min(partial):
                diff_before = cga[0]
            if min(cga[0]) > max(partial):
                diff_after = cga[0]
        elif cga_len == 2:  # this is the expected case
            diff_before = cga[0]
            diff_after = cga[1]
        else:
            print("Fatal error")
            exit(1)

        ret_dict = {
            "intersect_size": len(intersect),
            "diff_before_size": len(diff_before),
            "diff_after_size": len(diff_after),
            "intersect_size_perc": round(len(intersect) / full_size * 100, 1),
            "diff_before_size_perc": round(len(diff_before) / full_size * 100, 1),
            "diff_after_size_perc": round(len(diff_after) / full_size * 100, 1)
        }

        return ret_dict

    def prefilter_candidates(self):
        print("=> Pre-filtering candidates")
        # self._drop_overlaps_of_cds_5_ends()
        self._drop_unwanted_classes()
        # self._drop_cross_overlapping_annotations()

    def _drop_overlaps_of_cds_5_ends(self) -> None:
        cds_ref_df = self.gff_df[
            (self.gff_df["type"] == "gene")
            & (self.gff_df["attributes"].str.contains("gene_biotype=protein_coding;"))
        ].copy()
        fwd_cds_df = cds_ref_df[cds_ref_df["strand"] == "+"].copy()
        rev_cds_df = cds_ref_df[cds_ref_df["strand"] == "-"].copy()
        fwd_cds_df["end"] = fwd_cds_df["start"] + 30
        rev_cds_df["start"] = rev_cds_df["end"] - 30
        cds_ref_df = pd.concat([fwd_cds_df, rev_cds_df])
        cds_ref_bed = pybed.BedTool.from_dataframe(cds_ref_df).sort()
        tbl_bed = pybed.BedTool.from_dataframe(self.anno_tbl_df).sort()
        self.anno_tbl_df = tbl_bed.subtract(
            cds_ref_bed, s=True, F=0.99, A=True
        ).to_dataframe(names=self.gff_columns)

    def _drop_unwanted_classes(self) -> None:
        # drop unwanted overlaps like tRNA or rRNA
        drop_ref_df = self.gff_df[
            (self.gff_df["type"] == "tRNA")
            | (self.gff_df["attributes"].str.contains("gene_biotype=tRNA;"))
            | (self.gff_df["type"] == "rRNA")
            | (self.gff_df["attributes"].str.contains("gene_biotype=rRNA;"))
        ].copy()
        drop_ref_bed = pybed.BedTool.from_dataframe(drop_ref_df)
        # anno_tbl_bed = anno_tbl_bed.intersect(drop_ref_bed, v=True, f=0.10, s=True)
        tbl_bed = pybed.BedTool.from_dataframe(self.anno_tbl_df).sort()
        self.anno_tbl_df = tbl_bed.subtract(
            drop_ref_bed, s=True, F=0.10, A=True
        ).to_dataframe(names=self.gff_columns)

    # def _drop_cross_overlapping_annotations(self):
    #     ref_df = self.gff_df[self.gff_df["type"] == "gene"].copy()
    #     ref_bed = pybed.BedTool.from_dataframe(ref_df).sort()
    #     tbl_bed = pybed.BedTool.from_dataframe(self.anno_tbl_df).sort()
    #     tbl_df = tbl_bed.intersect(ref_bed, s=True, C=True, wa=True).to_dataframe(
    #         names=self.gff_columns + ["count"]
    #     )
    #     cross_df = tbl_df[tbl_df["count"] > 1]
    #     # TODO

    def get_intergenic_flanks(self):
        ref_gff_pb = pybed.BedTool.from_dataframe(self.gff_df[self.gff_df["type"] == "gene"]).sort()
        intergenic_pb = pybed.BedTool.from_dataframe(self.classes[(self.classes["attributes"].str.contains("annotation_class=intergenic")) |
                                                                  (self.classes["attributes"].str.contains("annotation_class=ncRNA")) |
                                                                  (self.classes["attributes"].str.contains("antisense_to_ncRNA_region"))]).sort()
        d_columns = self.gff_columns + [f"downstream_flank_{c}" for c in self.gff_columns] + ["downstream_flank_distance"]
        u_columns = self.gff_columns + [f"upstream_flank_{c}" for c in self.gff_columns] + ["upstream_flank_distance"]
        d_genes_df = intergenic_pb.closest(ref_gff_pb, D="a", io=True, iu=True).to_dataframe(names=d_columns)
        u_genes_df = intergenic_pb.closest(ref_gff_pb, D="a", io=True, id=True).to_dataframe(names=u_columns)
        intergenic_df = pd.merge(left=d_genes_df, right=u_genes_df, on=self.gff_columns, how="outer")
        intergenic_df["upstream_flank_distance"] = intergenic_df["upstream_flank_distance"].abs()
        for i in intergenic_df.index:
            for direction in ["up", "down"]:
                if intergenic_df.at[i, f"{direction}stream_flank_attributes"] == ".":
                    intergenic_df[f"{direction}stream_flank_distance"] = np.nan
                    continue
                attr = Helpers.parse_attributes(intergenic_df.at[i, f"{direction}stream_flank_attributes"])
                intergenic_df.at[i, f"{direction}stream_flank_gene_type"] = attr["gene_biotype"] if "gene_biotype" in attr.keys() else ""
                intergenic_df.at[i, f"{direction}stream_flank_gene_name"] = attr["name"] if "name" in attr.keys() else ""
                intergenic_df.at[i, f"{direction}stream_flank_orientation"] = "same" \
                    if intergenic_df.at[i, f"{direction}stream_flank_strand"] == \
                                  intergenic_df.at[i, "strand"] else "opposite"
        drop_cols = []
        for col in self.gff_columns:
            drop_cols.append(f"upstream_flank_{col}")
            drop_cols.append(f"downstream_flank_{col}")

        intergenic_df.drop(columns=drop_cols, inplace=True)
        self.classes = pd.merge(left=self.classes, right=intergenic_df, on=self.gff_columns, how="left")
        self.classes.fillna("", inplace=True)
        self.classes = Helpers.warp_non_gff_columns(self.classes)

    def get_gff_sequences_features(self, is_rna=False, slice_size=40) -> None:
        print("=> Extracting RNA features")
        gff_df = self.classes.copy()
        sequence_col = f"{'RNA' if is_rna else 'DNA'}_sequence"
        with tempfile.NamedTemporaryFile(mode="w") as temp:
            temp.write(self.fasta.full_fasta_str)
            fasta = pybed.example_filename(temp.name)
            temp.flush()
            gff_pb = pybed.BedTool.from_dataframe(gff_df).sort()
            seqs_str = (
                open(gff_pb.sequence(fi=fasta, rna=is_rna, s=True, tab=True).seqfn)
                .read()
                .replace(":", "\t")
                .replace("(-)", "r")
                .replace("(+)", "f")
                .replace("-", "\t")
                .replace("r", "\t-")
                .replace("f", "\t+")
            )
            os.remove(f"{temp.name}.fai")
        seqs_df = pd.read_csv(
            io.StringIO(seqs_str),
            sep="\t",
            names=[
                "seqid",
                "start",
                "end",
                "strand",
                sequence_col,
            ],
        )
        seqs_df["start"] += 1
        gff_df = pd.merge(
            left=gff_df,
            right=seqs_df,
            on=["seqid", "start", "end", "strand"],
            how="left",
        )
        gff_df.fillna("_", inplace=True)
        gff_df["GC_content"] = gff_df["RNA_sequence"].map(lambda x: round(SeqUtils.GC(x), 2))
        slice_prefix = f"last_{slice_size}_nt_"
        gff_df[f"{slice_prefix}RNA_sequence"] = gff_df["RNA_sequence"].str[-slice_size:]
        gff_df = Helpers.explode_dict_yielding_func_into_columns(gff_df, "RNA_sequence", self.get_rna_structure_scores)
        gff_df = Helpers.explode_dict_yielding_func_into_columns(gff_df, f"{slice_prefix}RNA_sequence", self.get_rna_structure_scores, slice_prefix)
        gff_df.drop(columns=["RNA_sequence", f"{slice_prefix}RNA_sequence"], inplace=True)
        self.classes = Helpers.warp_non_gff_columns(gff_df)

    def _get_poly_u_score(self, seq_str):
        ret_dict = {}
        seq_list = list(seq_str)
        u_indices = np.array([i for i, a in enumerate(seq_list, 1) if a == "U"])
        u_indices.sort()
        if u_indices.size == 0:
            return ret_dict
        #u_indices_groups = [list(group) for group in consecutive_groups(u_indices)]
        max_interrupt = 3
        u_stretches = np.split(u_indices, np.where(np.diff(u_indices) >= max_interrupt)[0] + 1)
        u_stretches_content = [len(s) for s in u_stretches]
        u_stretches_quality = [round(len(s)/(max(s) - min(s) + 1), 2) for s in u_stretches]
        print(u_stretches)
        print(u_stretches_content)
        print(u_stretches_quality)


        return ret_dict

    def get_rna_structure_scores(self, seq_str) -> dict:
        ret_dict = {}
        # create fold_compound data structure (required for all subsequently applied  algorithms)
        fc = RNA.fold_compound(seq_str)
        # compute MFE and MFE structure
        (mfe_struct, mfe) = fc.mfe()
        ret_dict["MFE"] = mfe
        # rescale Boltzmann factors for partition function computation
        fc.exp_params_rescale(mfe)
        # compute partition function
        (pp, pf) = fc.pf()
        ret_dict["partition_function"] = pf
        # compute centroid structure
        (centroid_struct, dist) = fc.centroid()
        ret_dict["centroid_structure_distance"] = dist
        # compute free energy of centroid structure
        centroid_en = fc.eval_structure(centroid_struct)
        ret_dict["centroid_MFE"] = centroid_en
        # compute MEA structure
        (MEA_struct, MEA) = fc.MEA()
        ret_dict["MEA"] = MEA
        # compute free energy of MEA structure
        MEA_en = fc.eval_structure(MEA_struct)
        ret_dict["MEA_MFE"] = MEA_en
        ret_dict["ensemble_MFE_structure_frequency"] = fc.pr_structure(mfe_struct)
        ret_dict["ensemble_diversity"] = fc.mean_bp_distance()

        return ret_dict
