import os

directory = 'RNAplfold_results'
directory_list = []
if os.path.isdir(directory):
    for path, dirs, files in os.walk(directory):
        for f in files:
            if '_dp.ps' in f: continue
            part = f.split('_')
            directory_list.append(os.path.abspath(path) + '/' + f)

to_remove = '/mnt/c/Irfan/Personal/Masters/BS6204/Project/RNAplfold_results/'
results_id = {}
FLANK = 50
LENGTH = 28
window = 50
with open('UnpairedProbability.50.txt', 'w') as write_file:
    write_file.write(f'id,{",".join([str(i+1) for i in range(LENGTH)])}\n')

    for k, file in enumerate(directory_list):
        seq_id = file.replace('_lunp','').replace(to_remove,'')
        # reading_id[count-1] = seq_id
        open_file = open(file, 'r')
        lines = open_file.readlines()

        result_list = []
        for i in range(2,len(lines)):
            line = lines[i].strip('\n')
            if line.startswith('#'):
                continue

            col = line.split('\t')
            pos = int(col.pop(0))

            ## Pos index 1 - 26
            if pos < (FLANK+1) or pos >= (FLANK + 1 + LENGTH): continue
            score = float(col[window-1]) 
            result_list.append(score)

        assert len(result_list) == LENGTH, len(result_list) 
        write_file.write(f'{seq_id},{",".join([str(i) for i in result_list])}\n')