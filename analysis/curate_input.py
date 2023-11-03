
file = 'combined_output.removeOutliers.training.csv'

WINDOW = 10
FLANK = 50
LENGTH = 28

write_file = open('BS6204_input.10WINDOW.csv', 'w')

headers = [
	'Id',
	'Input',
	'Quartile',
	'guideScore'
]

write_file.write(f'{",".join(headers)}\n')

def reverseC(seq):
    RT = {'A':'T','C':'G','G':'C','T':'A', 'N':'N'}
    reverseComplement = ''
    for i in seq:
        nt = RT.get(i)
        reverseComplement += nt
    return reverseComplement[::-1]

with open(file, 'r') as f:
	next(f)
	for line in f:
		line = line.strip('\n')
		col = line.split(',')
		seq_id = col[0]
		seq = col[2]
		target_sequence = col[3]
		start = FLANK-WINDOW
		stop = FLANK+LENGTH+WINDOW

		input_seq = target_sequence[start:stop]
		assert (len(input_seq)) == LENGTH+WINDOW*2
		assert reverseC(seq) in input_seq
		assert reverseC(seq) == input_seq[WINDOW:WINDOW+LENGTH]

		quartile = col[8]
		guideScore = col[9]
		write_file.write(f'{seq_id},{input_seq},{quartile},{guideScore}\n')