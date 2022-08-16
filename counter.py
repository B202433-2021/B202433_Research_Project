#!/usr/local/bin/python3

# counts the number of successful mopacs and gasteigers 

with open("successful_mopacs_2030.txt", "r") as infile:
     mopacs = infile.read()
     mopac_list = mopacs.split("\n")
     mopac_list = mopac_list[:-1]

with open("successful_gast_4181.txt", "r") as infile:
     gast = infile.read()
     gast_list = gast.split("\n")
     gast_list = gast_list[:-1]

print(mopac_list)
print(gast_list)

count = 0
for i in mopac_list:
   if i in gast_list:
       count += 1

print(count)
