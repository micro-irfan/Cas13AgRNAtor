

import os
from scipy.stats import pearsonr, spearmanr

directory = 'RNAplfold_results'
directory_list = []
if os.path.isdir(directory):
    for path, dirs, files in os.walk(directory):
        for f in files:
            if '_dp.ps' in f: continue
            part = f.split('_')
            directory_list.append(os.path.abspath(path) + '/' + f)

file = 'combined_output.removeOutliers.training.csv'
results_dict = {}
with open(file, 'r') as f:
    next(f)
    for line in f:
        col = line.strip('\n').split(',')
        seq_id = col[0]
        score = float(col[-1])
        quartile = int(col[-2].replace('Q',''))
        results_dict[seq_id] = (score, quartile)


LENGTH = 28
FLANK = 25
TOTAL = FLANK + LENGTH + FLANK
WINDOW = 50

to_remove = '/mnt/c/Irfan/Personal/Masters/BS6204/Project/RNAplfold_results/'
reading_id = {}
forward_matrix = [[[0 for _ in range(len(results_dict))] for _ in range(WINDOW)] for k in range(TOTAL)]
count = 0
for k, file in enumerate(directory_list):
    
    seq_id = file.replace('_lunp','').replace(to_remove,'')
    if seq_id not in results_dict.keys(): continue
    count += 1
    
    reading_id[count-1] = seq_id
    open_file = open(file, 'r')
    lines = open_file.readlines()

    for i in range(2,len(lines)):
        line = lines[i].strip('\n')
        if line.startswith('#'):
            continue

        col = line.split('\t')
        pos = int(col.pop(0))

        ## Pos index 1 - 26
        if pos < (FLANK+1): continue
        
        p = pos-FLANK-1
        for w in range(0, len(col)):
            offset = w // 2 
            # print (pos, p, offset, p-offset, w, count-1, col[w], file, len(col))
            if p-offset < 0: break
            if p-offset >= TOTAL: continue
            try:
                forward_matrix[p-offset][w][count-1] = float(col[w]) 
            except: 
                print (pos, p, offset, p-offset, w, count-1, col[w], file)
                print (len(results_dict))

print ('Finished Creating Matrix')
print ("Creating Corr Matrix")


import numpy as np
def correlation_value_from_3D_matrix(matrix, success):
    correlation_matrix1 = [[0 for _ in range(WINDOW)] for _ in range(TOTAL)]
    print (len(matrix), len(matrix[0]), len(matrix[0][0]))
    for i in range(len(matrix)): # Position
        for j in range(len(matrix[i])): # Window
            temp_list = matrix[i][j]
            temp_list = [float(i) for i in temp_list]
            correlation_matrix1[i][j] = pearsonr(success, temp_list)[0]
            pearsonr_high_corr  = abs(correlation_matrix1[i][j]) > 0.15 
            if pearsonr_high_corr:
                print ('Position:',i,'Window:',j,correlation_matrix1[i][j])
    return correlation_matrix1

success = []
quartile = []
for k,v in reading_id.items():
    
    success.append(results_dict[v][0])
    quartile.append(results_dict[v][1])

mat = correlation_value_from_3D_matrix(forward_matrix, success)


write_file_pearson  = open('pearsonr.RNAplfold.csv','w')
write_file_pearson.write('header,{}\n'.format(",".join(str(item) for item in [num for num in range(1, WINDOW+1)])))

for row, j in enumerate(mat):
    write_file_pearson.write('row_{},{}\n'.format(str(row+1),",".join(str(item) for item in j)))

'''
Logic Implemented
# if pos == 25:
#     if w < 2:
#         forward_matrix[p][w][k] = float(c) 

# if pos == 26:
#     if w < 2:
#         forward_matrix[p][w][k] = float(c) 
#     if w>=2 and w<4:
#         forward_matrix[p-1][w][k] = float(c) 

# if pos == 27:
#     if w < 2:
#         forward_matrix[p][w][k] = float(c) 
#     if w>=2 and w<4:
#         forward_matrix[p-1][w][k] = float(c) 
#     if w>=2 and w<4:
#         forward_matrix[p-2][w][k] = float(c) 

'''