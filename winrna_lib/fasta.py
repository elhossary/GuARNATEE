import glob
import os
from more_itertools import consecutive_groups
from Bio import SeqIO


class Fasta:
    def __init__(self, fasta_paths: list):
        self.fasta_paths = fasta_paths
        self.fwd_seqs = {}
        self.rev_seqs = {}
        self.organisms = {}
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
                specie_name = f"{seq_desc[1]}_{seq_desc[2]}"
                if specie_name in self.organisms.keys():
                    self.organisms[specie_name].append(f"{seq_desc[0]}")
                else:
                    self.organisms[specie_name] = [seq_desc[0]]
                self.fwd_seqs[seq_record.id] = str(seq_record.seq)
                self.rev_seqs[seq_record.id] = str(seq_record.reverse_complement().seq)[
                    ::-1
                ]  # 3` to 5` direction of the reverse complement
    """
    @staticmethod
    def count_bases(base: str, input_seq: str):
        return input_seq.count(base)

    @staticmethod
    def get_longest_consecutive_bases(base: str, input_seq: str):
        input_seq_lst = list(input_seq)
    """