import sys

import random

genePairs = []

with open(sys.argv[1],'r') as fr:
    for line in fr:
        genePairs.append(line.strip())
        
        
with open(sys.argv[2],'w') as fw:
    for i in random.sample(genePairs,int(sys.argv[3])):
        fw.write("%s\n"%(i))