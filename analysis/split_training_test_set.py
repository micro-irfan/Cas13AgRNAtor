import random

filename = 'combined_output.raw.csv'
training = 'combined_output.raw.training.csv'
test = 'combined_output.raw.test.csv'

training_file = open(training, 'w')
test_file = open(test, 'w')

training_count = 0
test_count = 0

output_count = 555
eighty_percent = int(0.8 * 555)
twenty_percent = 555 - eighty_percent
print ("Training Length:", eighty_percent) # 444
print ("Test Length:", twenty_percent) # 111

output_index = [i for i in range(1,output_count+1)]

test_list = random.choices(output_index, k=twenty_percent)

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

test_file.write(f'{",".join(headers)}\n')
training_file.write(f'{",".join(headers)}\n')

with open(filename, 'r') as f:
	next(f)
	for c, line in enumerate(f):
		if c in test_list:
			test_file.write(line)
		else:
			training_file.write(line)
