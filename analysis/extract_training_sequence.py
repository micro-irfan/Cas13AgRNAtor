

file = 'training.id.txt'

with open(file, 'r') as f:
    list_id = {i.strip('/n'):None for i in f}

combined = 'combined_output.raw.75.csv'
with open(combined, 'r') as f:
    next(f)
    for line in f:
        col = line.strip('\n').split(',')
        grna = col[0]
        seq = col[3]
        list_id[grna] = seq

with open('target_sequence.75.csv', 'w') as write_file:
    write_file.write('id,seq\n')
    for k,v in list_id.items():
        if not v: continue
        write_file.write(f'{k},{v}\n')