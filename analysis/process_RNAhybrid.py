import subprocess
from scipy.stats import pearsonr, spearmanr
import pingouin as pg
import pandas as pd

direct_repeat = 'GATTTAGACTACCCCAAAAACGAAGGGGACTAAAAC'

# Define the reverse_complement function to compute the reverse complement of a sequence
def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(sequence))

def get_RNAhyb_MFE(seq, add=False):
    if add:
        g = direct_repeat + seq
    else:
        g = seq

    cmd = f'RNAhybrid -c -s 3utr_human {reverse_complement(seq)} {g}'
    output = subprocess.check_output(cmd, shell=True, text=True)
    lines = output.strip().split(':')
    mfe = float(lines[4])
    return mfe

from collections import namedtuple

store = namedtuple('store', 'seqID seq quartile guideScore')

data = []
file = 'combined_output.removeOutliers.training.csv'
with open(file, 'r') as f:
    next(f)
    for line in f:
        line = line.strip('\n')
        col = line.split(',')
        seqID = col[0]
        seq = col[2]
        quartile = int(col[8].replace('Q',''))
        score = float(col[9])
        tmp = store(seqID, seq, quartile, score)
        data.append(tmp)

mfe_file = 'RNAfold.results.csv'
crRNA_mfe_list = []
with open(mfe_file, 'r') as f:
    next(f) 
    for line in f:
        line = line.strip('\n')
        col = line.split('\t')
        crRNA_mfe_list.append(float(col[3]))

write_file_1 = open('raw.RNAhybrid.2.csv', 'w')
write_file_2 = open('corr.RNAhybrid.2.txt', 'w')
write_file_2.write('start,stop,PCorr,SCorrS,ScorrQ,PartialCorr,pval\n')

LENGTH = 28
success = [d.guideScore for d in data]
quartile = [d.quartile for d in data]
add = False
for i in range(LENGTH-1):
    if i != 0: add = False
    for j in range(i+1, LENGTH):
        if j-i == 1: continue
        print (f'Analyzing Range Bp {i}-{j}')
        tmp = []
        
        for d in data:
            tmp_seq = d.seq[i:j]
            mfe = get_RNAhyb_MFE(tmp_seq, add)
            tmp.append((tmp_seq, mfe))

        string1 = ','.join([i[0] for i in tmp])
        string2 = ','.join([str(i[1]) for i in tmp])
        write_file_1.write(f'{i}-{j} Seq,{string1}\n')
        write_file_1.write(f'{i}-{j} Mfe,{string2}\n')

        mfe_result = [i[1] for i in tmp]

        df = pd.DataFrame({'Score': success, 'mfe': mfe_result, 'crRNAmfe': crRNA_mfe_list})
        partial_corr_result = pg.partial_corr(data=df, x='mfe', y='Score', covar='crRNAmfe')

        pearsonr_correlation_mfe  = pearsonr(success, mfe_result)[0]
        spearmanr_correlation_mfe1  = spearmanr(success, mfe_result)[0]
        spearmanr_correlation_mfe2  = spearmanr(quartile, mfe_result)[0]
        print (f'Correlation Range Bp {i}-{j}: {pearsonr_correlation_mfe}, {spearmanr_correlation_mfe1}, {spearmanr_correlation_mfe2}, {partial_corr_result["r"].values[0]}')
        write_file_2.write(f'{i},{j},{pearsonr_correlation_mfe},{spearmanr_correlation_mfe1},{spearmanr_correlation_mfe2},{partial_corr_result["r"].values[0]},{partial_corr_result["p-val"].values[0]}\n')



