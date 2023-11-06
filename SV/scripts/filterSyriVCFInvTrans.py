import re
import sys
pattern01 = "[A-Z]+"
pattern02 = "ChrB=[A-z]+[0-9]+"

pattern03 = "END=[0-9]+"

regexp01 = re.compile(pattern01)
regexp02 = re.compile(pattern02)
regexp03 = re.compile(pattern03)

sample_name = sys.argv[1]

add_line = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'

#print("OK")

fw = open(sys.argv[3],'w')
fr = open(sys.argv[2],'r')
fw02 = open(sys.argv[4],'w')

#i = 0
#print("OK")

for line in fr:
    if line.startswith("##fileformat="):
        fw.write("%s\n"%("##fileformat=VCFv4.1"))
    elif line.startswith("##"):
        fw.write(line)
    elif line.startswith("#CHROM"):
        fw.write("%s\n"%(add_line))
        temp_list = line.strip().split("\t")
        temp_list.append("FORMAT")
        temp_list.append(sample_name)
        fw.write("%s\n"%("\t".join(temp_list)))
    else:
        temp_list = line.strip().split("\t")
        if regexp01.findall(temp_list[2])[0] in ["INV","TRANS"] and regexp02.findall(temp_list[7])[0].split("=")[1] == temp_list[0]:
            #if (len(temp_list[3]) >= 50 and temp_list[3].count("N") == 0) or (len(temp_list[4]) >= 50 and temp_list[4].count("N") == 0):
            temp_list.append("GT")
            temp_list.append("1/1")
            fw.write("%s\n"%("\t".join(temp_list)))
            #i += 1
            
            var_type = regexp01.findall(temp_list[2])[0]
            var_len = int(regexp03.findall(temp_list[7])[0].replace("END=","")) - int(temp_list[1]) + 1
            fw02.write("%s\t%s\t%d\n"%(sample_name,var_type,var_len))
            

fr.close()
fw.close()
fw02.close()