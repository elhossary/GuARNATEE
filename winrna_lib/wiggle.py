import os
import re
import logging as logger
import pandas as pd
import numpy as np
from tqdm import tqdm
from Bio import Entrez


class Wiggle:
    def __init__(self, file_path, chrom_sizes=None):
        np.set_printoptions(suppress=True)
        self.file_path = file_path
        self.track_type = ""
        self.track_name = ""
        self.seqids = []
        self.coverages = {}
        self.signals = {}
        self.orientations = {}
        self.chrom_sizes = chrom_sizes if chrom_sizes is not None else {}
        self.spans = {}
        self.parse()
        self.generate_1d_signal()

    def parse(self):
        logger.info(f"Parsing {self.file_path}")
        current_wiggle_meta = {}
        with open(self.file_path, "r") as raw_file:
            logger.info(f"==> Loading file: {os.path.basename(self.file_path)}")
            file_header, all_contents, empty_seqids = self._parse_wiggle_str(
                raw_file.read()
            )
            current_wiggle_meta = self.parse_wiggle_header(
                file_header, current_wiggle_meta
            )
            for content_header, content in tqdm(
                all_contents.items(),
                desc=f"=> Loading wiggle file: {os.path.basename(self.file_path)}", bar_format='{desc} |{bar:20}| {percentage:3.0f}%'
            ):
                current_wiggle_meta = self.parse_wiggle_header(
                    content_header, current_wiggle_meta
                )
                seqid = current_wiggle_meta["variableStep_chrom"]
                chrom_size = (
                    self.chrom_sizes[seqid]
                    if seqid in self.chrom_sizes.keys()
                    else self.get_seq_length(seqid)
                )
                self.coverages[seqid] = np.fromstring(content, sep=" ").reshape(-1, 2)
                if chrom_size is None:
                    chrom_size = self.coverages[seqid][-1, 0]
                    logger.warning(
                        f"\nCould not get the chromosome size online, "
                        f"maximum wiggle position ({chrom_size}) will be used"
                    )
                else:
                    logger.info(f"====> Size retrieved for {seqid} is: {chrom_size}")
                self.chrom_sizes[seqid] = chrom_size
                self.spans[seqid] = current_wiggle_meta["variableStep_span"]
                self.orientations[seqid] = (
                    "-" if np.any(self.coverages[seqid][:, 1] < 0) else "+"
                )
            s = "\n     └── "
            logger.info(
                f"===> Parsed condition: {self.track_name}\n"
                f"     + Sequence IDs included:\n"
                f"     └── {s.join(self.coverages.keys())}"
            )

            self.seqids = list(self.coverages.keys()) + empty_seqids
            self.seqids.sort()
            if len(empty_seqids) > 0:
                logger.info(
                    f"     + Sequence IDs not included:\n"
                    f"     └── {s.join(empty_seqids)}"
                )

    def generate_1d_signal(self):
        for k in tqdm(self.coverages.keys(), desc="==> Parsing wiggles:", bar_format='{desc} |{bar:20}| {percentage:3.0f}%'):
            df1 = pd.DataFrame(np.arange(1, self.chrom_sizes[k] + 1), columns=["loc"])
            df2 = pd.DataFrame(self.coverages[k], columns=["loc", "score"])
            df2["loc"] = df2["loc"].astype(int)
            self.signals[k] = np.abs(
                pd.merge(df1, df2, right_on="loc", left_on="loc", how="outer")
                .fillna(0.0)["score"]
                .to_numpy()
            )

    @staticmethod
    def _parse_wiggle_str(in_str):
        ret_dict = {}
        header_text = in_str.split("\n", maxsplit=1)[0]
        in_str = in_str.replace(f"{header_text}\n", "")
        all_headers = re.findall(
            r"^.*chrom=.*$", in_str, flags=re.MULTILINE | re.IGNORECASE
        )
        splitters = ""
        for header in all_headers:
            splitters += header + "|"
        splitters = f"({splitters[:-1]})"
        split_str_list = re.split(
            rf"{splitters}", in_str, flags=re.MULTILINE | re.IGNORECASE
        )
        content_list = [i for i in split_str_list if i != ""]
        empty_seqids = []
        for i in range(0, len(content_list), 2):
            if content_list[i + 1].isspace():
                empty_seqids.append(
                    content_list[i].split("chrom=")[-1].split(" ")[0].replace('"', "")
                )
                continue
            ret_dict[content_list[i]] = content_list[i + 1]
        return header_text, ret_dict, empty_seqids

    def parse_wiggle_header(self, line, current_wiggle_meta):
        if "type=" in line:
            self.track_type = line.split("type=")[-1].split(" ")[0].replace('"', "")
        if "name=" in line:
            self.track_name = line.split("name=")[-1].replace("\n", "").replace('"', "")
        if "chrom=" in line:
            current_wiggle_meta["variableStep_chrom"] = (
                line.split("chrom=")[-1].split(" ")[0].replace('"', "")
            )
        if "span=" in line:
            current_wiggle_meta["variableStep_span"] = (
                line.split("span=")[-1].replace("\n", "").replace('"', "")
            )
        return current_wiggle_meta

    @staticmethod
    def get_seq_length(seq_id) -> int or None:
        Entrez.email = "muhammad_elhossary@hotmail.com"
        try:
            logger.info(f"===> Trying to get sequence ID length for: {seq_id}")
            esummary_handle = Entrez.esummary(
                db="nuccore", id=seq_id, report="full", idtype="acc"
            )
            esummary_record = Entrez.read(esummary_handle)[0]
        except Exception as e:
            logger.warning(f"===> Error fetching metadata: {e.args}")
            return None
        seq_len = int(
            str(esummary_record["Length"])
            .replace("IntegerElement(", "")
            .split(",", maxsplit=1)[0]
        )
        return seq_len
