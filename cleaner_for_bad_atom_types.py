#!/usr/local/bin/python3

# filters the complexes to remove ones where the pdbqt atom types were not consistent across docking 

with open("successful_both.txt", "r") as infile:
     successes = infile.read()
     success_list = successes.split("\n")
     success_list = success_list[:-1]

with open("wronguns.txt", "r") as infile:
     wronguns = infile.read()
     wrongun_list = wronguns.split("\n")
     wrongun_list = wrongun_list[:-1]

print(success_list)
print(wrongun_list)

with open("successful_both_clean.txt", "w") as outfile:
    count = 0
    for i in success_list:
        if i not in wrongun_list:
            count += 1
            outfile.write(i + "\n")
print(count)
