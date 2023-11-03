
file = 'combined_output.removeOutliers.training.csv'

list_id = []
with open(file, 'r') as f:
    next(f)
    for line in f:
        col = line.strip('\n').split(',')
        list_id.append(col[0])

with open('training.id.txt', 'w') as write_file:
    for i in list_id:
        write_file.write(f'{i}\n')