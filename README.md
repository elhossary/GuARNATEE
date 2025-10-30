# GuaRNAtee
## What is GuaRNAtee?
GuaRNAtee is a genome-wide computational tool for the suggestive identification of full-length small RNA transcripts in prokaryotes. It integrates differential RNA-Seq (dRNA-Seq) and Term-Seq coverage signals to detect coherent 5′ and 3′ transcript extremities for each condition. It relies on a cutoff-free, outlier-based approach applied in sliding windows over transcriptome coverage to extract extremities based solely on the local signal behavior. These paired coverage anomalies are nominated as potential sRNA loci, enabling the robust detection of transcript boundaries even under low signal-to-noise conditions that are typical for bacterial sRNAs. GuaRNAtee therefore prioritizes candidates for downstream validation and genome annotation without overclaiming final certainty.
In benchmarking against existing annotations, GuaRNAtee was able to recall almost all previously annotated small RNAs, while additionally suggesting numerous novel candidates with rankings. By directly leveraging sequencing-derived terminus evidence, GuaRNAtee supports hypothesis generation in small RNA biology and facilitates the discovery of previously unannotated RNA elements, including those embedded within coding sequences or UTRs.
### Disclaimer
This software is still a beta version

## Usage
### Prerequisites:
- Python >= 3.6
- Please refer to [requirement.txt](https://github.com/foerstner-lab/GuARNATEE/blob/main/requirements.txt) for the dependencies' installation.
### How to run:
- GuaRNAtee takes coverage files in wiggle format (.wig), reference genome in Fasta format, and reference annotations in GFF3 format.
- The ```--wigs``` argumant takes this format: ```<path>:<strand>:<condition>:<replicate>:<protocol>```, the fields are:
  - **path** path to WIG file
  - **strand** f or r
  - **condition** experimental growth condition (e.g. LB_0.4 M63)
  - **replicate** numeric replicate ID
  - **protocol** one of: ```TEX_pos```, ```TEX_neg```, or ```term```
- Use it as follows:
```lua
python3 guarnatee.py \
  --wigs <wig:strand:cond:rep:protocol> [more ...] \
  --gffs <gff> [more ...] \
  --fastas <fasta> [more ...] \
  --out_dir <output_dir>

```
- Example:
```shell
python3 guarnatee.py \
  --wigs \
    LB_0.4_rep1_forward.wig:f:LB_0.4:1:TEX_neg \
    LB_0.4_rep1_reverse.wig:r:LB_0.4:1:TEX_neg \
    LB_0.4_rep2_forward.wig:f:LB_0.4:2:TEX_pos \
    LB_0.4_rep2_reverse.wig:r:LB_0.4:2:TEX_pos \
    Term_LB_0.4_rep1_forward.wig:f:LB_0.4:1:term \
    Term_LB_0.4_rep1_reverse.wig:r:LB_0.4:1:term \
  --gffs ref_anno.gff known_ncRNAs.gff \
  --fastas refseq.fa \
  --out_dir results

```
## Future developments
Currently, GuaRNAtee is being further developed to detect sRNAs based on the coverage of conventional unfragmented sRNA-Seq libraries, as well as further improving the outlier detection algorithm.

## License
This software is licensed under [MIT license](https://github.com/foerstner-lab/GuARNATEE/blob/main/LICENSE)

## Citation
If you use this software, please refer to [CITATION.cff](https://github.com/foerstner-lab/GuaRNAtee/blob/main/CITATION.cff)