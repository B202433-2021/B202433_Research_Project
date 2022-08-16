#removes hydrogen from N-terminus of trimmed peptides
import sys
import numpy as np
arg= str(sys.argv[-1])
pdb_nonH=[]
pdb_H=[]
pdb_new=[]
num=[]
file=open('%s' %(arg),'r')
for x in file.readlines():
    if (x.split()[0]=='ATOM' and x.split()[2]!='H') or x.split()[0]=='HETATM':
        pdb_nonH.append(x)
    elif x.split()[0]=='ATOM' and x.split()[2]=='H':
        pdb_H.append(x)
file.close()
remove=[]
for x in pdb_nonH:
    #fix for resi more than 3 digits
    if len(x.split())==11:
      if x.split()[10]=='N1+':
        print (x)
        for y in range(len(pdb_H)):
            if pdb_H[y].split()[4]==x.split()[4]:
                remove.append(y)
                break
                
    else:
      if x.split()[11]=='N1+':
        print (x)
        for y in range(len(pdb_H)):
            if pdb_H[y].split()[5]==x.split()[5]:
                remove.append(y)
                break
file=open('%s' %(arg),'w')
for x in pdb_nonH:
    file.write(x)
for x in range(len(pdb_H)):
    if any(n==x for n in remove):
        continue
    else:
        file.write(pdb_H[x]) 
file.close()
