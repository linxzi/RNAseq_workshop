#! /usr/bin/env python3
import sys
import os

name = sys.argv[1]
matching_files = [f for f in os.listdir(".") if name in f]

#print(matching_files)

results = []
for input_file in matching_files:
    with open(input_file, 'r') as f:
        for line in f:
            if "#" not in line:
                results.append(line.strip().split())
                break

#print(results)

resultsCut = [[sublist[2], sublist[3], sublist[10]] for sublist in results]

#print(resultsCut)
n = 0
for filename in matching_files:
    last_underscore = filename.rfind("_")
    first_dot = filename.find(".", last_underscore)
    filename = filename[last_underscore + 1:first_dot]
    resultsCut[n].insert(0, filename)
    n += 1

#print(resultsCut)

resultsCut.insert(0, ["mat", "id", "alen", "eval"])

print(resultsCut)

data = resultsCut
# Print header
print(f"{data[0][0]:<10} {data[0][1]:<15} {data[0][2]:<20} {data[0][3]:<25}")
print("-" * 52)

# Print data rows
for row in data[1:]:
    print(f"{row[0]:<10} {row[1]:<15} {row[2]:<20} {data[0][3]:<25}")






