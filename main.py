from SmithWaterman import FilterSW_new


import os
from minineedle import needle, smith, core
from Bio import SeqIO
path = 'C:/Users/betya/BIOINF/miRNA_30_03_2023_summ_basecalled/pass'
t_file = 'fastq_runid_85e01dc3204a3a0a5e8e2cffab6a59b13e938a42_7_0.fastq'
all_files = os.listdir(path)
rest_files = len(all_files)
for file in all_files: 
    FilterSW_new(file)
    rest_files -= 1
    print(f"{file} был отсортирован по индексам, осталось {rest_files} файлов")