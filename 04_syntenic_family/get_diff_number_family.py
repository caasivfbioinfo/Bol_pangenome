#get gene family existing in 27,26,25... genomes 
# python 1.py numbers result_of_mSynOrths out

import sys

out=open(sys.argv[3],'w')
list1=[]
with open(sys.argv[2],'r') as fi:
    f = fi.readlines()
    for lines in f:
        list1=[]
        list2=[]
        line = lines.strip().split('\t')
        for i in line:
            b=i.split(':')[0]
            list2.append(i[:-1])
            if b not in list1:
                list1.append(b)
        if len(list1) == int(sys.argv[1]):
            out.write('\t'.join(list2) + '\n')