import sys
import re
from Bio import SeqIO

input_vcf = sys.argv[1]
input_ref = sys.argv[2]
output_vcf = sys.argv[3]

ref_dict = {}

for rec in SeqIO.parse(input_ref,'fasta'):
    ref_dict[rec.id] = str(rec.seq)


pattern = "[A-Z]+"
regexp = re.compile(pattern)

fr = open(input_vcf,'r')
fw = open(output_vcf,'w')

i = 0
j = 0

for line in fr:
    if line.startswith("##"):
        fw.write(line)
    elif line.startswith("#CHROM"):
        temp_fr = open(sys.argv[4],'r')
        for new_line in temp_fr:
            chr_id = new_line.split()[0]
            chr_len = new_line.split()[1]
            fw.write("##contig=<ID=%s,length=%s>\n"%(chr_id,chr_len))
            
        fw.write(line)
    elif not line.startswith("#"):
        var_type = regexp.findall(line.strip().split("\t")[2])[0]
        if var_type == "INS":
            chr_id = line.strip().split("\t")[0]
            location = int(line.strip().split("\t")[1])
            ref_base = line.strip().split("\t")[3]
            
            if ref_base == ref_dict[chr_id][location-1]:
                fw.write(line)
                i +=1
                
        elif var_type == "DEL":
            chr_id = line.strip().split("\t")[0]
            location = int(line.strip().split("\t")[1])
            ref_base = line.strip().split("\t")[3]
            
            if ref_base == ref_dict[chr_id][location-1:location-1+len(ref_base)]:
                fw.write(line)
                j += 1
                
print("Contained Ins: ",i)
print("Contained Del: ",j)


fr.close()
fw.close()