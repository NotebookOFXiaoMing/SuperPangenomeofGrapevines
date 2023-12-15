import sys
import itertools

fr = open(sys.argv[1],'r')

fw = open(sys.argv[2],'w')



total_id = []

for line in fr:
    id_list = line.strip().split()
    for i in list(itertools.combinations(id_list,2)):
        fw.write("%s\t%s\n"%(i[0],i[1]))
        # if i[0] not in total_id:
        #     total_id.append(i[0])
        # if i[1] not in total_id:
        #     total_id.append(i[1])
        
fw.close()


# fw = open(sys.argv[3],'w')
# fw.write("%s\n"%('\n'.join(total_id)))