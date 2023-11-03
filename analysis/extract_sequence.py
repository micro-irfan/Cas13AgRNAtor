
import numpy as np

def readfq(fp, gzipped = False): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                l = l.decode() if gzipped else l
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last =  last[1:], [], None #last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            l = l.decode() if gzipped else l
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield (name, ''.join(seqs), None) # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                l = l.decode() if gzipped else l
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield (name, seq, ''.join(seqs)) # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield (name, seq, None) # yield a fasta record instead
                break

def reverseC(seq):
    RT = {'A':'T','C':'G','G':'C','T':'A', 'N':'N'}
    reverseComplement = ''
    for i in seq:
        nt = RT.get(i)
        reverseComplement += nt
    return reverseComplement[::-1]

fasta = 'sequence.fasta'
gene = ['KRAS', 'PPIB', 'MALAT1', 'CLUC', 'GLUC']
sequence = {}

with open(fasta, 'r') as f:
    for name, seq, _ in readfq(f):
        for g in gene:
            if g in name:
                name = g
                break
        else:
            continue
        sequence[name] = seq

spacer = []

with open('sequence_list.txt', 'r') as f:
    for line in f:
        spacer.append(line.strip('\n').upper())

LENGTH = 28
FLANK  = 75

combined = 'combined_output.raw.75.csv'
combined = open(combined, 'w')

headers = [
    'Id', 
    'Gene',
    'Spacer Sequence',
    'Target Sequence', 
    'Pos',
    'Score1',
    'Score2',
    'Score3',
]

combined.write(f'{",".join(headers)}\n')

curated_spacer = set()

for g in gene:

    with open(f'{g.lower()}.scores.tsv', 'r') as file:
        cur_seq = sequence[g]

        next(file)
        for c, line in enumerate(file):
            if line.startswith('#'): continue

            line = line.strip('\n')
            col  = line.split('\t')

            index = int(col[0])
            score1 = col[1]
            score2 = col[2]
            score3 = col[3]

            diff_start, diff_end = 0, 0
            if g in ['GLUC', 'CLUC']:
                if index-LENGTH < 0:
                    index += 1


                start = index-LENGTH
                stop  = index
            else:
                index -= 1
                start = index
                stop = index+LENGTH

            print (index, 1, start, stop, stop-start, len(cur_seq))
            seq = cur_seq[start:stop]
            seq = reverseC(seq)

            assert len(seq) == LENGTH, (c, g, index, len(seq))
            assert seq in spacer, (c, g, index, seq, reverseC(seq))

            flank_start = start - FLANK
            
            if flank_start < 0:
                diff_start = 0 - flank_start
                flank_start = 0
            
            flank_end = stop + FLANK 
            if flank_end > len(cur_seq):
                diff_end = flank_end - len(cur_seq)
                flank_end = len(cur_seq)

            print (index, 2, flank_start, flank_end, diff_start, diff_end, flank_end - flank_start)
            flanking_sequence = cur_seq[flank_start:flank_end]
            flanking_sequence = 'N' * diff_start + flanking_sequence + "N" * diff_end

            combined.write(f'{g}-{c},{g},{seq},{flanking_sequence},{index-1},{score1},{score2},{score3}\n') 
            curated_spacer.add(seq)

            assert reverseC(seq) == flanking_sequence[FLANK:FLANK+LENGTH]
            assert len(flanking_sequence) == (FLANK + LENGTH + FLANK)

combined.close()

#missing 

missing = set(spacer) - curated_spacer
print (missing)
