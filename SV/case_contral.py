#python 1.py sv.vcf 1.list 2.list out
import sys

out = open(sys.argv[4],'w')

list1=[]
with open(sys.argv[2],'r') as fi:
    f = fi.readlines()
    for lines in f:
        line = lines.strip()
        list1.append(line)
list2=[]
with open(sys.argv[3],'r') as fi:
    f = fi.readlines()
    for lines in f:
        line = lines.strip()
        list2.append(line)
sam1=[]
sam2=[]
n =0
r = 0
with open(sys.argv[1],'r') as fi:
    f = fi.readlines()
    for lines in f:
        if lines[:2] != '##':
            line = lines.strip().split('\t')
            if line[0] == '#CHROM':
                for i in line[9:]:
                    n += 1
                    if i in list1:
                        sam1.append(n)
                    elif i in list2:
                        sam2.append(n)
            else:
                if ',' not in line[3] and ',' not in line[4]:
                    s = 0
                    CC1 = 0
                    CG1 = 0
                    GG1 = 0
                    NN1 = 0
                    CC2 = 0
                    CG2 = 0
                    GG2 = 0
                    NN2 = 0
                    for j in line[9:]:
                        s +=1
                        if s in sam1:
                            if j[:3] != './.':
                                r1 = j.split(':')[2].split(',')[0]
                                r2 = j.split(':')[2].split(',')[1]
                                if r1 != '.' and r2 != '.':
                                    if int(r1) >= 3 and int(r2) < 3:
                                        CC1 += 1
                                    elif int(r1) >= 3 and int(r2) >= 3:
                                        CG1 += 1
                                    elif int(r1) < 3 and int(r2) >= 3:
                                        GG1 += 1
                                    else:
                                        NN1 += 1
                            else:
                                NN1 += 1
                                    
                        if s in sam2:
                            if j[:3] != './.':
                                r1 = j.split(':')[2].split(',')[0]
                                r2 = j.split(':')[2].split(',')[1]
                                if r1 != '.' and r2 != '.':
                                    if int(r1) >= 3 and int(r2) < 3:
                                        CC2 += 1
                                    elif int(r1) >= 3 and int(r2) >= 3:
                                        CG2 += 1
                                    elif int(r1) < 3 and int(r2) >= 3:
                                        GG2 += 1
                                    else:
                                        NN2 += 1
                            else:
                                NN2 += 1
                        
                    out.write(line[2] + '\t' + str(CG1+GG1) + '\t' + str(CC1)+ '\t' + str(CG2+ GG2) + '\t' + str(CC2) + '\n')