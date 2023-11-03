import subprocess
import re
from scipy.stats import pearsonr, spearmanr

direct_repeat = 'GATTTAGACTACCCCAAAAACGAAGGGGACTAAAAC'
spacer_struct = '..(((((....((((.........)))).)))))..'

close_bracket_pos1 = set([c for c,i in enumerate(spacer_struct) if i =='('])
close_bracket_pos2 = set([c for c,i in enumerate(spacer_struct) if i ==')'])

write_file = open('RNAfold.results.csv', 'w')
write_file.write('Id\tquartile\tguideScore\tmfe\tgq\tdr\tresult\n')

def get_mfe(seqID, gRNA, quartile, guidescore):
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
    
    write_file.write(f'{seqID}\t{quartile}\t{guidescore}\t{mfe}\t{gq}\t{dr}\t{lines[1]}\n')
    print (seqID, guidescore, mfe, gq, dr, lines[1])
    return mfe, gq, dr

# Define the eval_fold function as needed
def eval_fold(output_line):
    # close_bracket_pos_results1 = set([c for c,i in enumerate(output_line) if i =='('])
    # close_bracket_pos_results2 = set([c for c,i in enumerate(output_line) if i ==')'])
    # result1 = not (close_bracket_pos1 - close_bracket_pos_results1)
    # result2 = not (close_bracket_pos2 - close_bracket_pos_results2)

    # if result1 and result2:
    #     return 1

    if output_line[:len(spacer_struct)] == spacer_struct:
        return 1
    
    return 0

from collections import namedtuple


store = namedtuple('store', 'mfe gq dr gRNA quartile score')
filename = 'combined_output.removeOutliers.training.csv'

results = []
with open(filename, 'r') as f:
    next(f)
    for line in f:
        line = line.strip('\n')
        col = line.split(',')

        seqID = col[0]
        gRNA = col[2]
        guidescore = float(col[9])
        quartile = int(col[8].replace('Q',''))

        mfe, gq, dr = get_mfe(seqID, gRNA, quartile, guidescore)
        tmp = store(mfe, gq, dr, gRNA, quartile, guidescore)
        results.append(tmp)

success = [v.score for v in results]
success1 = [v.quartile for v in results]

## Mfe        
mfe = [v.mfe for v in results]
gq = [v.gq for v in results]
dr = [v.dr for v in results]

print (gq)

spearmanr_correlation_mfe = spearmanr(success, mfe)[0]
pearsonr_correlation_mfe  = pearsonr(success, mfe)[0]

spearmanr_correlation_gq1  = spearmanr(success, gq)[0]
spearmanr_correlation_dr1  = spearmanr(success, dr)[0]

spearmanr_correlation_gq2  = spearmanr(success1, gq)[0]
spearmanr_correlation_dr2  = spearmanr(success1, dr)[0]

print (spearmanr_correlation_mfe, pearsonr_correlation_mfe, spearmanr_correlation_gq1, spearmanr_correlation_dr1)
print (spearmanr_correlation_mfe, pearsonr_correlation_mfe, spearmanr_correlation_gq2, spearmanr_correlation_dr2)

write_file.close()

# Example usage:
# result = get_mfe("sequence")
# print(result)

'''
0.1627139273798263 0.16316940593317703 -0.06552389769077079 -0.12831762954155565
0.1627139273798263 0.16316940593317703 -0.05075824753042268 -0.12474908888330456
'''