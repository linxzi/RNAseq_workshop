#! /usr/bin/env python3

delimiter = "^"
listCut = []
listSam = []
dictCounts = {}

with open("bowtie2.sam", "r") as file:
    

    for line in file:
        listSam.append(line.strip().split("\t"))

    #print(f"Reading SAM file...{listSam[0:5]}")
    i = 0
    for i in range(len(listSam)):
        listSam[i][2] = listSam[i][2].split(delimiter,1)[0]
        i += 1

    listCut = [[sublist[0], sublist[2]] for sublist in listSam]
    
    # n = 0
    # for line in listCut:
    #     n += 1
    #     print(f"{line}")
    #     if n == 30:
    #         break

    unique_col0_per_col1 = {}

    for row in listCut:
        col0_value = row[0]
        col1_value = row[1]

        if col1_value not in unique_col0_per_col1:
            unique_col0_per_col1[col1_value] = set()
        unique_col0_per_col1[col1_value].add(col0_value)

    result = {}
    for col1_value, unique_col0_values in unique_col0_per_col1.items():
        result[col1_value] = len(unique_col0_values)
    
    n = 0
    for key, value in result.items():
        n += 1
        print(f"{key}\t{value}")
        if n == 20:
            break







