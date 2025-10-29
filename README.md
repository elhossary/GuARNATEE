# GuARNATEE - Genome-wide Anomaly-based sRNA Transcripts Extremities Extraction
## What is GuARNATEE?
GuARNATEE is a genome-wide computational tool for the suggestive identification of full-length small RNA transcripts in prokaryotes. It integrates differential RNA-Seq (dRNA-Seq) and Term-Seq coverage signals to detect coherent 5′ and 3′ transcript extremities for each condition. It relies on a cutoff-free, outlier-based approach applied in sliding windows over transcriptome coverage to extract extremities based solely on the local signal behavior.

These paired coverage anomalies are nominated as potential sRNA loci, enabling the robust detection of transcript boundaries even under low signal-to-noise conditions that are typical for bacterial sRNAs. GuARNATEE therefore prioritizes candidates for downstream validation and genome annotation without overclaiming final certainty.

In benchmarking against existing annotations, GuARNATEE was able to recall almost all previously annotated small RNAs, while additionally suggesting numerous novel candidates. By directly leveraging sequencing-derived terminus evidence, GuARNATEE supports hypothesis generation in small RNA biology and facilitates the discovery of previously unannotated RNA elements, including those embedded within coding sequences or UTRs.
### Disclaimer
This software is still a beta version

## Usage

## Future developments
Currenly, GuARNATEE is being further developed to detect sRNAs based on the coverage conventional unfragmented sRNA-Seq libraries, and further improving the outlier detection algorithm

## License
This software is licensed under MIT license

## Citation
If you use this software, please refer to [CITATION.cff](https://github.com/elhossary/GuARNATEE)