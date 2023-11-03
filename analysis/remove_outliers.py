import numpy as np

## Create Dataset 

gene = ['KRAS', 'PPIB', 'MALAT1', 'CLUC', 'GLUC']
score_list = {g:[] for g in gene}

def find_closest_pair(numbers):
    if len(numbers) < 2:
        return []

    # Sort the NumPy array in ascending order
    numbers = np.sort(numbers)

    # Initialize variables to keep track of the closest pair and the minimum difference
    closest_pair = [numbers[0], numbers[1]]
    min_difference = abs(numbers[1] - numbers[0])

    # Iterate through the array to find the closest pair
    for i in range(len(numbers) - 1):
        difference = abs(numbers[i + 1] - numbers[i])
        if difference < min_difference:
            min_difference = difference
            closest_pair = [numbers[i], numbers[i + 1]]

    return closest_pair

def remove_outliers(id, data):
    data = np.array(data)

    # Calculate the mean and standard deviation
    mean = np.mean(data)
    std_dev = np.std(data)

    # Calculate Coefficient of Variance
    cv = (std_dev / mean)
    if cv > 0.1:
        print (id, 1, data, mean, std_dev, cv)

    # Set a threshold for Z-scores
    threshold = 1 # You can adjust this threshold as needed

    # Calculate Z-scores for the data points
    z_scores = (data - mean) / std_dev
    outlier = [abs(i) > threshold for i in z_scores]
    print (id, 2, data, z_scores, outlier)

    if any(outlier) and cv > 0.1:
        new_data = [data for outlier, data in zip(outlier, data) if not outlier]
        
        if len(new_data) == 1:
            print ("USING CLOSEST PAIR")
            new_data = find_closest_pair(data)
        
        # Calculate the mean and standard deviation
        mean = np.mean(new_data)
        std_dev = np.std(new_data)

        # Calculate Coefficient of Variance
        cv = (std_dev / mean)
        print ('ADJUSTED')
        print (id, 3, new_data, mean, std_dev, cv)
        if cv > 0.1:
            print ('OH NO')

        data = new_data

    return np.array(data)

#Id,Gene,Spacer Sequence,Target Sequence,Pos,Score1,Score2,Score3

file = 'combined_output.raw.training.csv'
open_file = 'adjusted_raw_data.training.csv'
write_file = open(open_file, 'w')
write_file.write('id,gene,mean,cv,std,pos,score1,score2,score3\n')
with open(file, 'r') as f:
    next(f)
    for line in f:
        line = line.strip('\n')
        col = line.split(',')
        g = col[1]

        scores = []
        replicate = 2
        scores.append(float(col[5]))
        scores.append(float(col[6]))
        if col[7]:
            scores.append(float(col[7]))
            replicate = 3

        scores = remove_outliers(col[0], scores)

        seq_id = col[0]
        tmp = list(scores)
        #write_file.write(f'{seq_id},{g},{",".join([str(i) for i in tmp])}\n')
        
        mean = np.mean(scores)
        std_dev = np.std(scores)
        cv = (std_dev / mean) * 100

        add_comma = "," if len(scores) == 2 else ""
        write_file.write(f'{seq_id},{g},{mean},{cv},{std_dev},{int(col[4])+2},{",".join([str(i) for i in tmp])}{add_comma}\n')

        mean = -np.log2(mean)
        score_list[g].append(mean)

        

# Define quartile values for the pooled dataset (P)
LQP = np.percentile(np.concatenate([v for k,v in score_list.items()]), 25)
MQP = np.percentile(np.concatenate([v for k,v in score_list.items()]), 50)
UQP = np.percentile(np.concatenate([v for k,v in score_list.items()]), 75)

print ('LQ pooled:', LQP)
print ('MQ pooled:', MQP)
print ('UQ pooled:', UQP)

lower_q = {g:0 for g in gene}
upper_q = {g:0 for g in gene}
median_q = {g:0 for g in gene}
for g,v in score_list.items():
    l_q = np.percentile(v, 25)
    m_q = np.percentile(v, 50)
    u_q = np.percentile(v, 75)

    print (f'Gene {g} LQ: {l_q}')
    print (f'Gene {g} MQ: {m_q}')
    print (f'Gene {g} UQ: {u_q}')

    upper_q[g] = u_q
    lower_q[g] = l_q
    median_q[g] = m_q

headers = [
    'Id', 
    'Gene',
    'Spacer Sequence',
    'Target Sequence', 
    'Pos',
    'unadjustedQuartile',
    'knockdownScore',
    'unadjustedlog2FC',
    'quartile',
    'adjustedlog2FC'
]

# Function to assign quartiles for each value
def assign_quartile(value, Q1, Q2, Q3):
    if value < Q1:
        return "Q1"
    elif Q1 <= value < Q2:
        return "Q2"
    elif Q2 <= value < Q3:
        return "Q3"
    else:
        return "Q4"

# Function to calculate normalized values
def calculate_normalized_value(x, LQD, UQD, LQP, UQP):
    return ((x - LQD) / (UQD - LQD)) * (UQP - LQP) + LQP

with open('combined_output.removeOutliers.training.csv', 'w') as write_file:
    write_file.write(f'{",".join(headers)}\n')
    
    with open(file, 'r') as f:
        next(f)
        for line in f:
            line = line.strip('\n')
            col = line.split(',')
            g = col[1]

            scores = []
            replicate = 2
            scores.append(float(col[5]))
            scores.append(float(col[6]))
            if col[7]:
                scores.append(float(col[7]))
                replicate = 3

            scores = remove_outliers(col[0], scores)
        
            knockdown_average = np.mean(scores)
            average = -np.log2(knockdown_average)
            # score_list[g].append()

            LQD = lower_q[g]
            UQD = upper_q[g]
            MQD = median_q[g]

            quartile = assign_quartile(average, LQD, MQD, UQD)

            adjusted_average = calculate_normalized_value(average, LQD, UQD, LQP, UQP)
            # if adjusted_average < 0: 
            #     adjusted_average = 0


            adjusted_quartile = assign_quartile(adjusted_average, LQP, MQP, UQP)

            # print (average, adjusted_average)
            write_file.write(f'{col[0]},{g},{col[2]},{col[3]},{col[4]},{quartile},{knockdown_average},{average},{adjusted_quartile},{adjusted_average}\n') 
