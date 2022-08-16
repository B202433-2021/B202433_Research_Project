#!/usr/local/bin/python3

# finds complexes where the pdbqt atom types were not consistent across the crystallographic, vina, gasteiger and mopac poses.

def get_complex_list():

    with open("atom_types.txt", "r") as infile:
        whole_file = infile.read()
        complexes = whole_file.split(">")[1:]
        return complexes

def dif_types(r,v,g,m):
#    print(r,v,g,m)
    for type in r: 
#        print(type)
        if r.count(type) != v.count(type):
            return True
        elif r.count(type) != g.count(type):
            return True
        elif r.count(type) != m.count(type):
            return True
    else:
        return False 
   

big_list = get_complex_list()

outfile = open("wronguns.txt", "w")
for comp in big_list:
#    print(comp)
    all_list = comp.split('\n')
#    print(all_list)
    pdb = all_list[0]
    ref = all_list[1]
    vina = all_list[2]
    gast = all_list[3]
    mopac = all_list[4]

#   print(f"Ref: {ref} \nVina: {vina} \nGast: {gast} \nMopac{mopac}")

    if dif_types(ref, vina, gast, mopac) :
        outfile.write(pdb + "\n")

outfile.close()
