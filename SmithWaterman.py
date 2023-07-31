from Bio import SeqIO
from minineedle import smith
from standart_primers import primers
path = r'C:\Users\betya\BIOINF\miRNA_30_03_2023_summ_basecalled\pass'
res_path = r'C:\Users\betya\BIOINF\miRNA_30_03_2023_summ_basecalled\results'

def FilterSW_new(filename):
    record_list = list(SeqIO.parse(fr"{path}/{filename}","fastq"))
    seq_list = [record_list[s].seq for s in range(len(record_list))]
    for no_seq in range(len(seq_list)):
        score = []
        for pr in range(1,21):
            aligment = smith.SmithWaterman((primers[f'primer{pr}']),seq_list[no_seq])
            score.append(aligment.get_score())
            max_score = max(score)
            max_index = score.index(max_score)
        if max_score >= 17:
            #matches[max_index].append(record_list[no_seq])
            with open(f"{path}/primers/primer{max_index}.fastq", "a") as output_handle:
                SeqIO.write(record_list[no_seq], output_handle, "fastq")
            
            
            
    return 




#record_list[1].format("fastq")
#record_list = list(SeqIO.parse(fr"{path}/{filename}","fastq"))
#id_list = [record_list[k].id for k in range(len(record_dict))]
#description_list = [record_list[s].description for s in range(len(record_dict))]
#quality_list = [record_list[q].letter_annotations["phred_quality"] for q in range(len(record_dict))]
#record_dict = SeqIO.index(fr"{path}/{filename}","fastq")