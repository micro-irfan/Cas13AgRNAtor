import subprocess
import re

direct_repeat = 'GATTTAGACTACCCCAAAAACGAAGGGGACTAAAAC'
spacer_struct = '..(((((....((((.........)))).)))))..'

close_bracket_pos1 = set([c for c,i in enumerate(spacer_struct) if i =='('])
close_bracket_pos2 = set([c for c,i in enumerate(spacer_struct) if i ==')'])

write_file = open('RNAfold.results.all.csv', 'w')
write_file.write('Id\tquartile\tguideScore\tmfe\tgq\tdr\tresult\n')

# Define the eval_fold function as needed
def eval_fold(output_line):
    if output_line[:len(spacer_struct)] == spacer_struct:
        return 1
    
    return 0

def get_mfe(seqID, gRNA):
    crRNA = direct_repeat + gRNA
    # You may need to change the path to your RNAfold executable
    cmd = f'echo "{crRNA}" | RNAfold --gquad --noPS'
    output = subprocess.check_output(cmd, shell=True, text=True)
    
    # Parse the output to extract MFE, gq, and dr values
    lines = output.strip().split('\n')
    # print (lines)
    mfe_match = re.search(r'\(([^)]+)\)$', lines[1])
    mfe = float(mfe_match.group(1))
    gq = 1 if '+' in lines[1] else 0
    dr = eval_fold(lines[1])
    
    write_file.write(f'{seqID}\t{mfe}\t{gq}\t{dr}\t{lines[1]}\n')
    print (seqID, mfe, gq, dr, lines[1])
    return mfe, gq, dr

filename = 'combined_output.raw.75.csv'

results = []
with open(filename, 'r') as f:
    next(f)
    for line in f:
        line = line.strip('\n')
        col = line.split(',')

        seqID = col[0]
        gRNA = col[2]

        mfe, gq, dr = get_mfe(seqID, gRNA)