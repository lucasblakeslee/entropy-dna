import os
from pathlib import Path

def main():
    db_directory_names = list(Path(".").glob("*_databases"))
    for db_directory_name in db_directory_names:
        average_files(db_directory_name)

def average_files(directory_name):
    files = os.listdir(directory_name)
    sequences = {}
    for file_name in files:
        with open(os.path.join(directory_name, file_name)) as file_data:
            for line in file_data.readlines():
                sequence = line[:line.find(" ")]
                occurences = int(line[line.find(" ")+1:])
                if sequences.get(sequence):
                    sequences[sequence] = (sequences[sequence] + occurences)/2
                else:
                    sequences[sequence] = occurences

    with open(f"{directory_name}_averages.db_txt", "w+") as f:
        output = ""
        for key, value in sequences.items():
            output = output + f"{key} {value}\n"
        f.write(output[:-1])

if __name__ == "__main__":
    main()