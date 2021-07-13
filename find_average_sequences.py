import os

def average_files(directory_path):
    files = os.listdir(directory_path)

    sequences = {}
    for file_name in files:
        with open(directory_path + "/" + file_name) as file_data:
            for line in file_data.readlines():
                sequence = line[:line.find(" ")]
                occurences = int(line[line.find(" ")+1:])
                if sequences.get(sequence):
                    sequences[sequence] = (sequences[sequence] + occurences)/2
                else:
                    sequences[sequence] = occurences
    
    with open("gammaproteobacteria_averages.db_txt", "w+") as f:
        output = ""
        for key, value in sequences.items():
            output = output + f"{key} {value}\n"
        f.write(output[:-1])

average_files('gammaproteobacteria_databases')
