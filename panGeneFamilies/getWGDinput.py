import sys
import random
import re

pattern = "[0-9]+"
regexp = re.compile(pattern)

def myfun(input_list):
    new_list = []
    for i in input_list:
        if int(i) <= 1:
            new_list.append(0)
        else:
            new_list.append(1)
    if sum(new_list) == 0:
        return True
    else:
        return False


family_ids = []

with open(sys.argv[1],'r') as fr:
    for line in fr:
        family_ids.append(line.strip())
        
#family_ids = random.sample(family_ids,100)

singleCopyFamilyID = []

with open(sys.argv[2],'r') as fr:
    for line in fr:
        if line.startswith("OG") and line.split()[0] in family_ids and myfun(regexp.findall(line)[1:-1]):
            #print(regexp.findall(line)[1:-1])
            singleCopyFamilyID.append(line.split()[0])
                    
#print(singleCopyFamilyID)

fr = open(sys.argv[3],'r')
fw = open(sys.argv[4],'w')

for line in fr:
    if line.startswith("OG") and line.split()[0] in singleCopyFamilyID:
        species_list = line.strip().split("\t")[1:]
        gene_list = []
        for species in species_list:
            for i in species.split(", "):
                if len(i)>0 and i.count("_") == 0:
                    gene_list.append(i)

        
        #print("\t".join(output_list))
        if len(gene_list) >=2:
            print(len(gene_list))
            fw.write("%s\n"%("\t".join(gene_list)))
        
fw.close()