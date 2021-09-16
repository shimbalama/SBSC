# Description

Test resources

* `HCC1937_*_BRCA1.pileup.gz` are pileups from long-read-sequencing bams on
  HCC1937 tumour and control samples, selected in the BRCA1 region
  (chr17:43044294-43125363 in GRCh38 coords). The .tbi files are their indices.

* `chr17_BRCA1.fa.gz` is the reference fasta file for chr17, with all bases
  outside BRCA1 set to N so it still has the correct length but compresses down
  nicely to a few hundred kb.

* `out.tsv` is the expected variants called on the pileups.
