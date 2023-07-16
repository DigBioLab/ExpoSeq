import random
import string
import os

def generate_random_sequence(length):
    return ''.join(random.choice('ATCG') for _ in range(length))

def generate_random_quality(length):
    return ''.join(random.choice(string.ascii_letters) for _ in range(length))

def generate_fastq_file(file_path, num_records, sequence_length):
    with open(file_path, 'w', newline='\n') as file:
        for i in range(num_records):
            sequence = generate_random_sequence(sequence_length)
            quality = generate_random_quality(sequence_length)
            file.write(f"@Sequence{i+1}\n{sequence}\n+\n{quality}\n")

# Usage example:
generate_fastq_file('fastq3.fastq', 1000, 100)
