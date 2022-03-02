import glob
import io
import os
import pandas as pd
import pybedtools as pybed
from more_itertools import consecutive_groups
from Bio import SeqIO



class Fasta:
    def __init__(self, fasta_paths: list):
        self.fasta_paths = fasta_paths
        self.fwd_seqs = {}
        self.full_fasta_str = ""
        self.rev_seqs = {}
        self.organism_seqid_groups = {}
        self.parse()

    def parse(self):
        print("=> Parsing input fasta files")
        parsed_paths = []
        for item in self.fasta_paths:
            for sub_item in glob.glob(item):
                parsed_paths.append(os.path.abspath(sub_item))
        for fasta_path in parsed_paths:
            parse_tmp = SeqIO.parse(os.path.abspath(fasta_path), "fasta")
            for seq_record in parse_tmp:
                seq_desc = seq_record.description.split()
                specie_name = (
                    f"{seq_desc[1]}_{seq_desc[2]}"
                    if len(seq_desc) >= 3
                    else "undefined"
                )
                if specie_name in self.organism_seqid_groups.keys():
                    self.organism_seqid_groups[specie_name].append(f"{seq_desc[0]}")
                else:
                    self.organism_seqid_groups[specie_name] = [seq_desc[0]]
                self.fwd_seqs[seq_record.id] = str(seq_record.seq)

                # self.rev_seqs[seq_record.id] = str(seq_record.reverse_complement().seq)[
                #    ::-1
                # ]  # 3` to 5` direction of the reverse complement
        self.full_fasta_str = "\n".join(
            [f">{k}\n{v}" for k, v in self.fwd_seqs.items()]
        )

        print(f"==> Parsed {len(parsed_paths)} Fasta files")
