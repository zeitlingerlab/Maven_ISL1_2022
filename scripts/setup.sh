#!bin/bash

#Melanie Weilert
#July 2021
#Purpose: Obtain data and set up samples for running.
cd /n/projects/mw2098/collaboration/for_srivastava/data/bam/individual

wget https://www.dropbox.com/sh/05cxqaubheaxy96/AACnf3TYYYD_nddP8XTxycZDa/BWA_aligned_to_hg19/I_nko_d6cm_rep1_BWA_aligned.bam?dl=0
wget https://www.dropbox.com/sh/05cxqaubheaxy96/AABN7Wn-KO_OMcxn3o_nbd5La/BWA_aligned_to_hg19/I_nko_d6cm_rep2_BWA_aligned.bam?dl=0
wget https://www.dropbox.com/sh/05cxqaubheaxy96/AADRjD-fdJqtTNPlr1wer2NJa/BWA_aligned_to_hg19/I_nko_d6cm_rep3_BWA_aligned.bam?dl=0
wget https://www.dropbox.com/sh/05cxqaubheaxy96/AAAIdxCXfZnFyd380knXB2RTa/BWA_aligned_to_hg19/I_wt_d6cm_rep1_BWA_aligned.bam?dl=0
wget https://www.dropbox.com/sh/05cxqaubheaxy96/AAC-Vnko8vRU0qVluUVrHIJ9a/BWA_aligned_to_hg19/I_wt_d6cm_rep2_BWA_aligned.bam?dl=0
wget https://www.dropbox.com/sh/05cxqaubheaxy96/AAC36-5uK0e-0nSvXDpjFFpta/BWA_aligned_to_hg19/I_wt_d6cm_rep3_BWA_aligned.bam?dl=0
wget https://www.dropbox.com/sh/05cxqaubheaxy96/AAC7AXDAiyWIGHlTKmoE7kfMa/BWA_aligned_to_hg19/I_wt_d6cm_rep4_BWA_aligned.bam?dl=0
wget https://www.dropbox.com/sh/05cxqaubheaxy96/AAB8qMxaXIvuAeUjBxFc9x6Ta/BWA_aligned_to_hg19/I_wt_s3mn_rep1_BWA_aligned.bam?dl=0
wget https://www.dropbox.com/sh/05cxqaubheaxy96/AACZrI8QwbKzB7bboQ8MHOM9a/BWA_aligned_to_hg19/I_wt_s3mn_rep2_BWA_aligned.bam?dl=0
wget https://www.dropbox.com/sh/05cxqaubheaxy96/AAAROr9dPGfe8ixQceO0SuTWa/BWA_aligned_to_hg19/I_wt_s3mn_rep3_BWA_aligned.bam?dl=0
wget https://www.dropbox.com/sh/05cxqaubheaxy96/AAA-m1wpAHFVIU0t73XTPFnJa/BWA_aligned_to_hg19/input_nko_d6cm_rep1_BWA_aligned.bam?dl=0
wget https://www.dropbox.com/sh/05cxqaubheaxy96/AABbJ04CnNDuEFUzQ5it86qla/BWA_aligned_to_hg19/input_wt_d6cm_rep1_BWA_aligned.bam?dl=0
wget https://www.dropbox.com/sh/05cxqaubheaxy96/AACLMMqVM9tsCSz74S6RKuaha/BWA_aligned_to_hg19/input_wt_d6cm_rep2_BWA_aligned.bam?dl=0
wget https://www.dropbox.com/sh/05cxqaubheaxy96/AACgaR_RijR7j3DInXKGzHUwa/BWA_aligned_to_hg19/input_wt_d6cm_rep3_BWA_aligned.bam?dl=0
wget https://www.dropbox.com/sh/05cxqaubheaxy96/AAC0-6mXvdOOiR4zM1V8lke-a/BWA_aligned_to_hg19/input_wt_d6cm_rep4_BWA_aligned.bam?dl=0
wget https://www.dropbox.com/sh/05cxqaubheaxy96/AACRDE4mUhMYuDbnUsXjSy7Ga/BWA_aligned_to_hg19/input_wt_s3mn_rep1_BWA_aligned.bam?dl=0
wget https://www.dropbox.com/sh/05cxqaubheaxy96/AAC55jPPBBIRDgPLs1ZQG4yka/BWA_aligned_to_hg19/input_wt_s3mn_rep2_BWA_aligned.bam?dl=0
wget https://www.dropbox.com/sh/05cxqaubheaxy96/AADwYWE0OrpoNUlx_6ntTB4Ua/BWA_aligned_to_hg19/input_wt_s3mn_rep3_BWA_aligned.bam?dl=0
wget https://www.dropbox.com/sh/05cxqaubheaxy96/AAAK-ocfcTm0OGjBcwrqKAKaa/BWA_aligned_to_hg19/readme.txt?dl=0

wget https://www.dropbox.com/sh/05cxqaubheaxy96/AACwu6lJ8jjQBP35o2NsyDy2a/BWA_aligned_to_hg19/Other_BWA_aligned_to_hg19/input_iko_d6cm_rep1.bam?dl=0
wget https://www.dropbox.com/sh/05cxqaubheaxy96/AAAgtNoL3NzYA4k3xND3ZIiJa/BWA_aligned_to_hg19/Other_BWA_aligned_to_hg19/N_iko_d6cm_rep2.bam?dl=0
wget https://www.dropbox.com/sh/05cxqaubheaxy96/AADN3zguM6v3O-bHddEE3f2pa/BWA_aligned_to_hg19/Other_BWA_aligned_to_hg19/N_nko_d6cm_rep2.bam?dl=0
wget https://www.dropbox.com/sh/05cxqaubheaxy96/AADllRa2Hz-QU54i8GnQPvLCa/BWA_aligned_to_hg19/Other_BWA_aligned_to_hg19/N_nko_d6cm_rep3.bam?dl=0
wget https://www.dropbox.com/sh/05cxqaubheaxy96/AADh516JdTccZJXFJc044M56a/BWA_aligned_to_hg19/Other_BWA_aligned_to_hg19/N_wt_d6cm_rep1.bam?dl=0
wget https://www.dropbox.com/sh/05cxqaubheaxy96/AABKuaHLhZpl_MZQP-m0nk5-a/BWA_aligned_to_hg19/Other_BWA_aligned_to_hg19/N_wt_d6cm_rep2.bam?dl=0
wget https://www.dropbox.com/sh/05cxqaubheaxy96/AADLJj3DeqTht14_EMqyQg3_a/BWA_aligned_to_hg19/Other_BWA_aligned_to_hg19/N_wt_d6cm_rep3.bam?dl=0

#Rename to remove the dropbox extensions.
rename .bam?dl=0 .bam *.bam?dl=0
rename .txt?dl=0 .txt *.txt?dl=0
