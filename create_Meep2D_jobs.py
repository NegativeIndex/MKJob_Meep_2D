#!/usr/bin/env python3
import os
import re
import subprocess
import glob
import numpy as np
import shutil 
import sys

sys.path.insert(0,'/Users/wdai11/python-study/MeepFunctions')
import my_helper as my
import flatten
import itertools

######################### 
# global variables
##########################
class common:
    fname="detector.py"     
    otherfiles=()
    cpu=8
    memory=20

######################### 
# generate one job
##########################
def find_line_numbers(lines,regexp):
    numbers=[]
    for i,line in enumerate(lines):
        # print(line)
        if re.match(regexp, line):
            numbers.append(i)
    return numbers

######################### 
# generate one job
##########################
def other_steps(dname,cpu,memory):
    
    # step 3.1 copy other files
    for ff in common.otherfiles:
        print('Copy file '+ff)
        shutil.copy(ff,dname)

    # step 3.2, generate job file
    os.chdir(dname)
    fname=common.fname
    command=os.path.join(os.getenv("HOME"),
                         'python-study/MKJob_Meep_2D/mkjob-meep-2D')
    subprocess.call([command,fname])
   
    jobfiles=glob.glob('dwt*.job')
    assert len(jobfiles)==1,"One job file in a folder"
    with open(jobfiles[0], 'r') as f:
        lines=f.readlines()

  
    for i,line in enumerate(lines):
        line=re.sub(r'^#\$ -pe smp\s+[0-9]+',
                    '#$ -pe smp {}'.format(cpu), line)
        line=re.sub(r'^mpirun -np [0-9]+',
                  'mpirun -np {}'.format(cpu), line)
        line=re.sub(r'^#\$ -l mt=[0-9]+G',
                    '#$ -l mt={}G'.format(memory), line)
        lines[i]=line
        
    fout=open(jobfiles[0], 'wt')
    fout.write("".join(lines))
    fout.close()   
    # step 3.3 generate job.begin
    open("job.begin", 'a').close()
    print("Generate job.begin")

######################### 
# check two files
##########################
def check_two_files():
    with open(common.fname) as f:
        file1=f.read().splitlines()

    with open("create_Meep2D_jobs.py") as f:
        file2=f.read().splitlines()

    file1=[line.lstrip().rstrip() for line in file1 
           if re.search('=',line)]
    file2=[line.lstrip().rstrip()
           for line in file2 if re.search('=',line)]
    
    sameline=[line for line in file1 if line in file2]

    print('-----------------------')
    for line in sameline:
        print(line)
    print('-----------------------')
    userInput = input('The two files have such commone lines. y/n?  ');

    if userInput!='y':
        sys.exit('Something is wrong. Modify the file')
    
########################## 
#### main function
##########################
# collect job file parameters
check_two_files()

# prepare py file
path=os.path.abspath("./")
fname=common.fname

lines=[]
with open(fname, "rt") as fin:
    lines=fin.readlines()

line1=find_line_numbers(lines,'class common')
line2=find_line_numbers(lines,r'^##########################')

lmin=line1[0]
i=0
while line2[i]<lmin: i+=1
lmax=line2[i]
############################
# all the parameters
pols=["Ez",]
geoms=["Empty","RC"]
fcens=np.arange(0.245,0.31,0.03)
dfs=np.ones(fcens.size)*0.03

freqs=list(zip(fcens,dfs))
gparas=[(a,)
        for a in np.array([0.3,]) ]


paras1=list(itertools.product(pols, ("Empty",), freqs, 
                              ( (0,)*len(gparas[0]), )   ))
paras1=[flatten.flattenArrayN(elem) for elem in paras1 ]

paras2=list(itertools.product(pols, geoms[1:], freqs, gparas ))
paras2=[ flatten.flattenArrayN(elem) for elem in paras2 ]

paras=paras1+paras2
# paras=paras1

# [print(elem) for elem in paras]

######################
njobs=len(paras)
print('{} jobs are generated'.format(njobs))
userInput = input('Any key please, Ctrl-c to quit')

count=0
for para in paras:
    os.chdir(path)
    
    pol,geom,fcen,df,a=para
    sig='{}_{}_F{:5.3f}'.format(geom,pol,fcen)
    
    dname=sig
    print(sig)
    if not os.path.exists(dname):
        # step 1, build the folder
        os.makedirs(dname)
        # step 2, write the new py simulation file file
        fname2=os.path.join(dname,fname)
        for i in range(lmin,lmax):
            line=lines[i];
            # print(line)
            line=re.sub(r'^(\s+geom=).*$',r'\g<1>"{:s}"'.format(geom), line)
            line=re.sub(r'^(\s+pol=).*$',r'\g<1>"{:s}"'.format(pol), line)
            line=re.sub(r'^(\s+fcen=).*$', r'\g<1>{:0.3f}'.format(fcen), line)
            line=re.sub(r'^(\s+df=).*$', r'\g<1>{:0.3f}'.format(df), line)
            # line=re.sub(r'^(\s+a=).*$', r'\g<1>{:0.3f}'.format(a), line)
            # print(line)
            lines[i]=line

        fout=open(fname2, 'wt')
        fout.write("".join(lines))
        fout.close()
        print("Generate {:s}".format(fname2))

        # step 3, other steps
        other_steps(dname, common.cpu, common.memory)
        count+=1
        print('Generated {} jobs, {} to do'.format(count, njobs-count))
        print('-------------------')

print('{} jobs are generated'.format(count))

