import os
import io
import sys
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
    def __init__(self, gff_obj, anno_tbl_df: pd.DataFrame, fasta: Fasta):
        if anno_tbl_df.empty or gff_obj.gff_df.empty:
            print("Error: No candidates to classify")
            sys.exit(1)
        prohibited_types = ["region", "exon", "sequence_feature"]
        self.gff_df = gff_obj.gff_df[~gff_obj.gff_df["type"].isin(prohibited_types)].copy()
        self.anno_tbl_df = anno_tbl_df
        self.fasta = fasta
        self.gff_columns = gff_obj.column_names
        self.seqids = set.intersection(set(self.gff_df["seqid"].unique().tolist()),
                                       set(self.fasta.fwd_seqs.keys()),
                                       set(self.anno_tbl_df["seqid"].unique().tolist()))
        self.anno_tbl_df = self.anno_tbl_df[self.anno_tbl_df["seqid"].isin(self.seqids)]  # .iloc[:100, :]
        self.classes = pd.DataFrame()
        self.dispatch()  # run the classifier at init

    def dispatch(self):
        #self.prefilter_candidates()
        self.classify()
        self.get_intergenic_flanks()
        self.get_gff_sequences_features(is_rna=True)
        self._drop_redundancies()
        self.classes.sort_values(["seqid", "start", "end"], inplace=True)
        #print(self.classes.head(100).to_string())

    def _drop_redundancies(self):
        df = Helpers.expand_attributes_to_columns(self.classes)
        combs = list(
            product(df["seqid"].unique().tolist(), ["+", "-"], ["start", "end"])
        )
        drop_ids = []
        for seqid, strand, select_column in tqdm(
            combs, desc="===> Cleaning redundancies", bar_format='{desc} |{bar:20}| {percentage:3.0f}%'
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
        anno_tbl_df = Helpers.warp_non_gff_columns(self.anno_tbl_df)
        candidates_df = anno_tbl_df.copy()
        candidates_df["type"] = "sRNA_candidate"

        ref_df = Helpers.merge_same_intervals(self.gff_df)
        attr_filters = ["biotype", "name", "locus_tag", "old_locus_tag", "ID"]
        candidates_pb = pybed.BedTool.from_dataframe(
            Helpers.warp_non_gff_columns(candidates_df)
        ).sort()
        ref_pb = pybed.BedTool.from_dataframe(ref_df).sort()
        ref_columns = [f"ref_{x}" for x in self.gff_columns]
        intersect_columns = self.gff_columns + ref_columns + ["overlap_size"]
        candidates_df = candidates_pb.intersect(ref_pb, wao=True, f=0.5).to_dataframe(names=intersect_columns)

        candidates_df["overlap_size"] = candidates_df["overlap_size"].astype(int)
        candidates_df.drop(columns=["ref_seqid", "ref_source", "ref_score", "ref_phase"], inplace=True)
        ref_cols = [c for c in candidates_df.columns if "ref_" in c]
        # clean and prepare
        substr_filter_func = \
            lambda in_str, substrs, sep: \
                f"{sep}".join([w for w in in_str.split(sep) if any(substr in w for substr in substrs)])
        for rc in ref_columns:
            if rc not in candidates_df.columns:
                continue
            candidates_df.loc[candidates_df[rc] == ".", rc] = ""

        # get novel intergenic (no intersection or unannotated regions)
        diff_df = candidates_df[candidates_df["overlap_size"] <= 0].copy()
        candidates_df.drop(diff_df.index, inplace=True)
        diff_df.drop(columns=ref_cols + ["overlap_size"], inplace=True)
        diff_df["attributes"] += ";annotation_class=intergenic;detection_status=novel"
        self.classes = pd.concat([self.classes, diff_df], ignore_index=True)
        del diff_df
        #########################
        # Clean ref_type column
        candidates_df["ref_type_tmp"] = candidates_df["ref_type"].apply(substr_filter_func, args=(["RNA", "CDS"], "|")).replace("", np.nan)
        candidates_df.loc[candidates_df["ref_attributes"].str.contains("gene_biotype=protein_coding"), ["ref_type_tmp"]] = "CDS"
        candidates_df["ref_type"].update(candidates_df["ref_type_tmp"])
        candidates_df.drop(columns=["ref_type_tmp"], inplace=True)
        # filter attributes to get values of interest
        candidates_df = Helpers.filter_attributes(candidates_df, attr_filters, "ref_attributes")
        # Elemenate redundancy in candidates by grouping overlaps
        candidates_df = candidates_df.groupby(self.gff_columns, as_index=False).agg({c: list for c in ref_cols + ["overlap_size"]})

        ##############################
        # Split overlapping candidates into single/multiple overlaps
        single_intersect_df = candidates_df[candidates_df["ref_type"].map(len) == 1].copy()
        candidates_df.drop(single_intersect_df.index, inplace=True)
        multi_intersect_df = candidates_df.copy()
        del candidates_df
        for col in ref_cols + ["overlap_size"]:
            # unpack lists with single value while preserving data type
            single_intersect_df[col] = single_intersect_df[col].map(max)

        antisense_single_df = single_intersect_df.loc[single_intersect_df["strand"] != single_intersect_df["ref_strand"]].copy()
        single_intersect_df.drop(antisense_single_df.index, inplace=True)
        antisense_single_df["attributes"] = antisense_single_df["attributes"] \
                                            + ";annotation_class=antisense_to_" \
                                            + antisense_single_df["ref_type"] \
                                            + "_region;detection_status=novel;" \
                                            + antisense_single_df["ref_attributes"]
        antisense_single_df["overlap_fragment_ratio"] = \
            antisense_single_df["overlap_size"] / \
            (antisense_single_df["ref_end"] - antisense_single_df["ref_start"] + 1)\
            * 100
        antisense_single_df["overlap_fragment_ratio"] = antisense_single_df["overlap_fragment_ratio"].round(2)
        antisense_single_df.drop(columns=ref_cols + ["overlap_size"], inplace=True)
        antisense_single_df = Helpers.rewrap_attributes_column(antisense_single_df)
        antisense_single_df = Helpers.warp_non_gff_columns(antisense_single_df)
        self.classes = pd.concat([self.classes, antisense_single_df], ignore_index=True)
        del antisense_single_df
        ##########################################
        single_intersect_df = single_intersect_df.apply(self.add_overlap_info, args=[False], axis=1, result_type='expand')
        sorf_mask = (single_intersect_df["ref_type"] == "CDS") & (single_intersect_df["overlap_fragment_ratio"] > 75)
        orf_int_mask = (single_intersect_df["ref_type"] == "CDS") & (single_intersect_df["overlap_fragment_ratio"] <= 75)
        single_intersect_df.loc[sorf_mask, "attributes"] += ";annotation_class=sORF;detection_status=known"
        single_intersect_df.loc[orf_int_mask, "attributes"] += ";annotation_class=ORF_int;detection_status=novel"
        single_intersect_df.loc[sorf_mask, "attributes"] = \
            single_intersect_df.loc[sorf_mask, "attributes"] + ";" + single_intersect_df.loc[sorf_mask, "ref_attributes"]
        single_intersect_df.loc[orf_int_mask, "attributes"] = \
            single_intersect_df.loc[orf_int_mask, "attributes"] + ";" + single_intersect_df.loc[orf_int_mask, "ref_attributes"]
        tmp_df = single_intersect_df[single_intersect_df["ref_type"] != "CDS"].copy()
        tmp_df["attributes"] = tmp_df["attributes"]\
                               + ";annotation_class=" \
                               + tmp_df["ref_type"] \
                               + ";detection_status=known;" \
                               + tmp_df["ref_attributes"]
        single_intersect_df.update(tmp_df)
        del tmp_df
        single_intersect_df["start"] = single_intersect_df["start"].astype(int)
        single_intersect_df["end"] = single_intersect_df["end"].astype(int)
        single_intersect_df.drop(columns=ref_cols + ["overlap_size"], inplace=True)
        single_intersect_df = Helpers.rewrap_attributes_column(single_intersect_df)
        single_intersect_df = Helpers.warp_non_gff_columns(single_intersect_df)
        self.classes = pd.concat([self.classes, single_intersect_df], ignore_index=True)
        del single_intersect_df
        ######################################
        for i in multi_intersect_df.index:
            # clean cross double stranded intersections
            drop_ids = [x for x, y in enumerate(multi_intersect_df.at[i, "ref_strand"]) if y != multi_intersect_df.at[i, "strand"]]
            for ref_col in ref_cols + ["overlap_size"]:
                lst = multi_intersect_df.at[i, ref_col]
                multi_intersect_df.at[i, ref_col] = [y for x, y in enumerate(lst) if x not in drop_ids]
            if all(x == "CDS" for x in multi_intersect_df.at[i, "ref_type"]):
                multi_intersect_df.at[i, "ref_attributes"] = ";".join(multi_intersect_df.at[i, "ref_attributes"])
                multi_intersect_df.at[i, "attributes"] += f";annotation_class=cross_CDS;detection_status=novel;{multi_intersect_df.at[i, 'ref_attributes']}"
            elif any(x == "CDS" for x in multi_intersect_df.at[i, "ref_type"]) and any(x == "ncRNA" or "ORF_int" in x for x in multi_intersect_df.at[i, "ref_type"]):
                multi_intersect_df.at[i, "ref_attributes"] = ";".join(multi_intersect_df.at[i, "ref_attributes"])
                multi_intersect_df.at[i, "attributes"] += f";annotation_class=ORF_int;detection_status=known;{multi_intersect_df.at[i, 'ref_attributes']}"
            else:
                multi_intersect_df.at[i, "ref_attributes"] = ";".join(multi_intersect_df.at[i, "ref_attributes"])
                types = "_and_".join(multi_intersect_df.at[i, "ref_type"])
                multi_intersect_df.at[i, "attributes"] += f";annotation_class=cross_{types};detection_status=known;{multi_intersect_df.at[i, 'ref_attributes']}"
            for v_id, v in enumerate(multi_intersect_df.at[i, "overlap_size"]):
                ref_start = multi_intersect_df.at[i, "ref_start"][v_id]
                ref_end = multi_intersect_df.at[i, "ref_end"][v_id]
                ref_type = multi_intersect_df.at[i, "ref_type"][v_id]
                start = multi_intersect_df.at[i, "start"]
                end = multi_intersect_df.at[i, "end"]
                overlap_info = self._get_overlap_position((start, end), (ref_start, ref_end))
                rename_ref = "" if ref_type == "CDS" else f"{ref_type}_"
                multi_intersect_df.at[i, f"{rename_ref}overlap_fragment_ratio"] = overlap_info["intersect_size_perc"]
                if ref_type != "CDS":
                    continue
                if multi_intersect_df.at[i, "strand"] == "-":
                    overlap_info["diff_before_size_perc"], overlap_info["diff_after_size_perc"] = \
                        overlap_info["diff_after_size_perc"], overlap_info["diff_before_size_perc"]
                multi_intersect_df.at[i, f"{rename_ref}upstream_fragment_ratio"] = overlap_info["diff_before_size_perc"]
                multi_intersect_df.at[i, f"{rename_ref}downstream_fragment_ratio"] = overlap_info["diff_after_size_perc"]
        multi_intersect_df.drop(columns=ref_cols + ["overlap_size"], inplace=True)
        multi_intersect_df.fillna("", inplace=True)
        multi_intersect_df = Helpers.rewrap_attributes_column(multi_intersect_df)
        multi_intersect_df = Helpers.warp_non_gff_columns(multi_intersect_df)
        self.classes = pd.concat([self.classes, multi_intersect_df], ignore_index=True)
        del multi_intersect_df

    def add_overlap_info(self, row: pd.Series, intersection_only=False) -> pd.Series:
        if row["overlap_size"] == 0:
            return row
        overlap_info = self._get_overlap_position((row["start"], row["end"]), (row["ref_start"], row["ref_end"]))
        row["overlap_fragment_ratio"] = overlap_info["intersect_size_perc"]
        if intersection_only:
            return row
        if row["strand"] == "-":
            overlap_info["diff_before_size_perc"], overlap_info["diff_after_size_perc"] = \
                overlap_info["diff_after_size_perc"], overlap_info["diff_before_size_perc"]
        row["upstream_fragment_ratio"] = overlap_info["diff_before_size_perc"]
        row["downstream_fragment_ratio"] = overlap_info["diff_after_size_perc"]
        #row["gene_info"] = Helpers.get_naming_attributes(row["ref_attributes"], prefix=f"{row['ref_type']}_")
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
            sys.exit(1)

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
        intergenic_pb = pybed.BedTool.from_dataframe(self.classes[self.classes["attributes"].str.contains("annotation_class=intergenic")]).sort()
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

        slice_prefix = f"{slice_size}_nt_"
        gff_df[f"{slice_prefix}RNA_sequence"] = gff_df["RNA_sequence"].str[-slice_size:]
        gff_df[f"10nt_RNA_sequence"] = gff_df["RNA_sequence"].str[-10:]
        gff_df = Helpers.explode_dict_yielding_func_into_columns(gff_df, "RNA_sequence", self.get_rna_structure_scores)
        gff_df = Helpers.explode_dict_yielding_func_into_columns(gff_df, f"{slice_prefix}RNA_sequence", self.get_rna_structure_scores, slice_prefix)
        gff_df = Helpers.explode_dict_yielding_func_into_columns(gff_df, f"10nt_RNA_sequence", self._get_poly_u_score)
        gff_df.drop(columns=["RNA_sequence", f"{slice_prefix}RNA_sequence", "10nt_RNA_sequence"], inplace=True)
        self.classes = Helpers.warp_non_gff_columns(gff_df)

    def _get_poly_u_score(self, seq_str):
        ret_dict = {}
        seq_list = list(seq_str)
        u_indices = np.array([i for i, a in enumerate(seq_list, 1) if a == "U"])
        u_indices.sort()
        if u_indices.size == 0:
            ret_dict["poly_u_score"] = 0
            return ret_dict
        #u_indices_groups = [list(group) for group in consecutive_groups(u_indices)]
        """
        max_interrupt = 3
        u_stretches = np.split(u_indices, np.where(np.diff(u_indices) >= max_interrupt)[0] + 1)
        u_stretches_content = [len(s) for s in u_stretches]
        u_stretches_quality = [round(len(s)/(max(s) - min(s) + 1), 2) for s in u_stretches]
        """
        ret_dict["poly_u_score"] = len(u_indices) / len(seq_list)
        return ret_dict

    def get_rna_structure_scores(self, seq_str) -> dict:
        ret_dict = {}
        # create fold_compound data structure (required for all subsequently applied  algorithms)
        fc = RNA.fold_compound(seq_str)
        #min_pairs = 1
        #max_pairs = int(len(seq_str) / 2)
        #print(fc.(structure=fc))
        #print(fc.E_ext_hp_loop(min_pairs, max_pairs))
        # compute MFE and MFE structure
        (mfe_struct, mfe) = fc.mfe()
        base_pairing_nts = int(mfe_struct.count("(") + mfe_struct.count(")"))
        #ret_dict["base_pairs_count"] = base_pairing_nts
        ret_dict["base_pairs_content"] = round(base_pairing_nts / len(seq_str) * 100, 2)
        ret_dict["MFE"] = mfe
        """
        # rescale Boltzmann factors for partition function computation
        fc.exp_params_rescale(mfe)
        # compute partition function
        (pp, pf) = fc.pf()
        ret_dict["partition_func"] = pf
        # compute centroid structure
        (centroid_struct, dist) = fc.centroid()
        ret_dict["centroid_struct_dist"] = dist
        # compute free energy of centroid structure
        centroid_en = fc.eval_structure(centroid_struct)
        ret_dict["centroid_MFE"] = centroid_en
        # compute MEA structure
        (MEA_struct, MEA) = fc.MEA()
        ret_dict["MEA"] = MEA
        # compute free energy of MEA structure
        MEA_en = fc.eval_structure(MEA_struct)
        ret_dict["MEA_MFE"] = MEA_en
        ret_dict["ensemble_MFE_struct_freq"] = fc.pr_structure(mfe_struct)
        ret_dict["ensemble_diversity"] = fc.mean_bp_distance()
        """
        return ret_dict
