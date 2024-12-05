# Guaymas Loop Seq

Script use to process LoopSeq full length 16S rRNA gene sequence data for the manuscript in revision:
J.E. Hinkle, J.P. Chanton, M.A. Moynihan, S.E. Ruff, A. Teske, Complex bacterial diversity of Guaymas Basin hydrothermal sediments revealed by synthetic long-read sequencing (LoopSeq).

# Description

[Script](https://github.com/moyn413/GuaymasLoopSeq/blob/main/dada2_loopseq) used to process LoopSeq data using DADA2 and generate a [phyloseq](https://joey711.github.io/phyloseq/) object. This script follows the [DADA workflow](https://benjjneb.github.io/dada2/)) with modifications suggested for [long read sequencing](https://benjjneb.github.io/LRASManuscript/LRASms_fecal.html) and [Loop Genomics](https://github.com/benjjneb/dada2/issues/991). This workflow uses [cutadapt](https://cutadapt.readthedocs.io/en/stable/) to identify full-length 16S rRNA gene sequences. Reads without both primers are removed, as they are not considered full-length. However, note that it is possible to have relatively long sequence reads (e.g. 1200bp) that are not full-length and might still be useful for some analyses. As alternative to using cutadapt, sequences could be filtered by length.
