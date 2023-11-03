from scipy.stats import pearsonr, spearmanr
import itertools
import os

gRNA_LENGTH = 28
FLANK  = 50
LENGTH = FLANK + gRNA_LENGTH + FLANK

write_file = open('high_corr.RNAContext2.txt', 'w')

def count_func(seq, letter):
    c = 0
    length = 1
    for i in range(len(seq)-length):
        tmp = seq[i:i+length]
        assert len(tmp) == length
        if tmp.upper() in letter:
            c += 1
    return c

def hypothetical_max_count(seq_len, letter):
    if len(set(letter)) == 1:
        return seq_len - len(letter) + 1

    if letter == letter[::-1]:
        tmp = letter
        count = 1
        while len(tmp) < seq_len:
            count += 1
            tmp += letter[1:]

        if len(tmp) > seq_len:
            count -= 1

        return count
    
    return seq_len // len(letter)

def get_nucleotide_density(target_dict, success, quartile, analysis_length, letter = 'A'):
    #set_up
    mid_pos = int(FLANK/2)
    letter_length = 1
    #w = window length for sliding analysis

    forward  = [[[0 for i in range(LENGTH - FLANK + 1 - letter_length)] for j in range(FLANK)] for k in range(analysis_length)]
    print ('Processing {}'.format(letter))

    
    for i in range(mid_pos, LENGTH - mid_pos + 1 - letter_length):
        for j in range(0, mid_pos+1):
            for e in range(2):
                if e == 0 and j != 0:
                    tmp_len = (i+j+letter_length-1) - (i-j)
                else:
                    tmp_len = (i+j+letter_length) - (i-j)

                sequence_length = tmp_len # hypothetical_max_count(tmp_len, letter)
                for count, (k, seq) in enumerate(target_dict.items()):
                    if e == 0 and j != 0:
                        tmp = seq[i-j:i+j+letter_length-1]
                    else:
                        if j*2 == FLANK: continue
                        tmp = seq[i-j:i+j+letter_length]

                    # count_forward = count_func(tmp, letter)
                    count_forward = tmp.count(letter[0])
                    if len(letter) > 1:
                        count_forward += tmp.count(letter[1])
                    assert count_forward <= sequence_length
                    forward[count][len(tmp) - letter_length][i-mid_pos] = count_forward/sequence_length
                    #print (k, i, j, i-mid_pos, len(tmp)- letter_length, tmp, count_forward, sequence_length)

    file_name = f'Nucleotide_context_{letter}_edited.csv'
    pearson_corr, spearman_corr, keep = correlation_value_from_3D_matrix(forward, success, quartile, letter)
    if not keep:
        return

    write_file_spearman = open(f'./RNA_context2/spearmanr.{file_name}','w')
    write_file_pearson  = open(f'./RNA_context2/pearsonr.{file_name}','w')
    write_file_spearman.write('header,{}\n'.format(",".join(str(item) for item in [num for num in range(1, LENGTH-FLANK+1)])))
    write_file_pearson.write('header,{}\n'.format(",".join(str(item) for item in [num for num in range(1, LENGTH-FLANK+1)])))
    
    for row, j in enumerate(spearman_corr):
        write_file_spearman.write('row_{},{}\n'.format(str(row+1),",".join(str(item) for item in j)))
    for row, j in enumerate(pearson_corr):
        write_file_pearson.write('row_{},{}\n'.format(str(row+1),",".join(str(item) for item in j)))
    
    write_file_pearson.close()
    write_file_spearman.close()

def correlation_value_from_3D_matrix(matrix, success, quartile, letter):
    letter_length = 1
    pearsonr_correlation_matrix  = [[0 for i in range(LENGTH - FLANK + 1 - letter_length)] for i in range(FLANK)]
    spearmanr_correlation_matrix = [[0 for i in range(LENGTH - FLANK + 1 - letter_length)] for i in range(FLANK)]
    keep = True
    for i in range(len(matrix[0])): #50
        for j in range(len(matrix[0][i])): #>79
            temp_list = []
            for k in range(len(matrix)): 
                temp_list.append(matrix[k][i][j])

            spearmanr_correlation_matrix[i][j] = spearmanr(quartile, temp_list)[0]
            pearsonr_correlation_matrix[i][j]  = pearsonr(success, temp_list)[0] 

            if j > 25 and j <= 53:
                r = [(i,j) for i, j in zip(quartile, temp_list)]
                #print (letter, i, j, spearmanr_correlation_matrix[i][j])
                #print (r)

            spearmanr_high_corr = abs(spearmanr_correlation_matrix[i][j]) > 0.2 
            pearsonr_high_corr  = abs(pearsonr_correlation_matrix[i][j]) > 0.2 
            high_corr = spearmanr_high_corr or pearsonr_high_corr
            if high_corr:
                keep = True
                if abs(spearmanr_correlation_matrix[i][j]) > 0.2:
                    print ('spearmanr : {}, {}, range:{}, pos:{}'.format(letter, spearmanr_correlation_matrix[i][j], i+1, j+1)) 
                    write_file.write('spearmanr : {}, {}, range:{}, pos:{}'.format(letter, spearmanr_correlation_matrix[i][j], i+1, j+1) + '\n')
                if abs(pearsonr_correlation_matrix[i][j]) > 0.2:
                    print ('pearsonr  : {}, {}, range:{},  pos:{}'.format(letter, pearsonr_correlation_matrix[i][j], i+1, j+1)) 
                    write_file.write('pearsonr  : {}, {}, range:{},  pos:{}'.format(letter, pearsonr_correlation_matrix[i][j], i+1, j+1) + '\n')
                
    return pearsonr_correlation_matrix, spearmanr_correlation_matrix, keep

filename = 'combined_output.removeOutliers.training.csv'
x = ['A', 'C', 'G', 'T']
x = ['A', 'C', 'G', 'T', ('A','T'), ('C','G'), ('A','C'), ('G','T')] # , ('A','C'), ('G','T')
# perm_2 = [''.join(p) for p in itertools.product(x, repeat=2)]
# perm_3 = [''.join(p) for p in itertools.product(x, repeat=3)]
# perm_4 = [''.join(p) for p in itertools.product(x, repeat=4)]
# perm_5 = [''.join(p) for p in itertools.product(x, repeat=5)]
# perm_6 = [''.join(p) for p in itertools.product(x, repeat=6)]
new_list = x # + perm_2 + perm_3 + perm_4  # + perm_5 # + perm_6
target_dict = {}
success = []
quartile = []
with open(filename, 'r') as f:
    next(f)
    for line in f:
        line = line.strip('\n')
        col  = line.split(',')
        success.append(float(col[9]))
        quartile.append(int(col[8].replace('Q','')))
        target_dict[col[0]] = col[3]

analysis_length = len(target_dict)
print ("Finish Parsing throung RNA context file!")

# get_nucleotide_density(target_dict, success, analysis_length, 'AA')
for i in new_list:
  get_nucleotide_density(target_dict, success, quartile, analysis_length, i)

write_file.close()