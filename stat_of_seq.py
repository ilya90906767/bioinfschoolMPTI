import seaborn as sns
import os
from Bio import SeqIO
from minineedle import smith
from standart_primers import primers
path = 'C:/Users/betya/BIOINF/miRNA_30_03_2023_summ_basecalled/pass'

def FilterSW_new(filename):
    record_list = list(SeqIO.parse(fr"{path}/{filename}","fastq"))
    seq_list = [record_list[s].seq for s in range(len(record_list))]
    matches = [[] for _ in range(40)]
    for no_seq in range(len(seq_list)):
        score = []
        for pr in range(1,21):
            aligment = smith.SmithWaterman((primers[f'primer{pr}']),seq_list[no_seq])
            score.append(aligment.get_score())
            max_score = max(score)
            max_index = score.index(max_score)
        if max_score >= 17:
            matches[max_index].append(record_list[no_seq])
            
    quality_list = [record_list[q].letter_annotations["phred_quality"] for q in range(len(record_list))]        
    return matches,quality_list


test = FilterSW_new(os.listdir(path)[0])
sns.displot(test[1][11], kde=True, bins=15) # статистика одного выровненного рида 
max_list = []                                 #список номеров выровненных ридов, которые имеют качество от 70 до 95
for k in range(len(test[1])):
    max_list.append(max(test[1][k]))
for z in range(len(max_list)):
    if (max_list[z]<95) and (max_list[z]>70):
        print(z)
sns.displot(max_list, kde=True, bins=15)  #график максимумов q в выравнивании на один праймер