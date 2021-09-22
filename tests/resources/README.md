# Description

Test resources

* `HCC1937_*_RNF157.pileup.gz` are pileups from long-read-sequencing bams on
  HCC1937 tumour and control samples, selected in the RNF157 gene region
  (chr17:76142456-76240493 in GRCh38 coords). The .tbi files are their indices.

* `chr17_RNF157.fa.gz` is the reference fasta file for chr17, with all bases
  outside RNF157 set to N so it still has the correct length but compresses down
  nicely to a few hundred kb.

* `output_001_9_12_10.tsv` is the expected variants called on the pileups.
          ^   ^ ^  ^
          |   | |  min_depth_normal
          |   | min_depth_tumour
          |   min_mean_base_quality
          P_value 
