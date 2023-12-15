import sys


with open(sys.argv[1],'r') as fr:
    i = 0
    seq_seq = ""
    for line in fr:
        i += 1
        if i == 1:
            seq_ids = line.strip().split("-")
        else:
            seq_seq += line.strip()
            
    with open(sys.argv[2],'w') as fw:
        
        fw.write(">%s\n%s\n"%(seq_ids[0],seq_seq[0:int(len(seq_seq)/2)]))
        
        fw.write(">%s\n%s\n"%(seq_ids[1],seq_seq[int(len(seq_seq)/2):]))