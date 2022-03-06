
for i in range(141):
    with open(f'{i}.fasta', 'r') as file:
        line = file.readline()
        with open(f'{i}.fasta', 'w') as file:
            print('>', file=file)
            print(line, file=file)