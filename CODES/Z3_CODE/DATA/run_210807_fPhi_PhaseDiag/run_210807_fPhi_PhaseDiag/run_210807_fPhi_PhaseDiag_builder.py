#! /usr/bin/env python

import subprocess
import numpy
import os

#partition_info=['normal',56]; time_str='1-00:00:00' # = [partition,ncores]
partition_info=['small',56]; time_str='2-00:00:00' # = [partition,ncores]
#partition_info=['large',24] # = [partition,ncores]
#partition_info=['development', 56]; time_str='02:00:00' # = [partition,ncores]; time_str='2-00:00:00'
project_name=os.getcwd().split('/')[-3]
myemail=os.environ["MYEMAIL"]
charge_str="inverted_mbl"

#Log the current submission
logstr='''run_210807_fPhi_PhaseDiag:
For L=18: if l == 18 and G<0.3 and vnum%4==0 and Ee<0:



***write me***

Codebase: 
Cluster: ls5
gittag: '''

#Do the git-ing
cmd='git commit -a -m "Commit before run run_201108_1_juliakrylov" > /dev/null'
print(cmd)
subprocess.call(cmd,shell=True)

cmd='gittag.py > temp.temp'
subprocess.call(cmd,shell=True)
logstr+=open('temp.temp','r').readline()
subprocess.call('rm temp.temp',shell=True)

open('run_210807_fPhi_PhaseDiag.log','w').write(logstr)

#Setup the versionmap and qsub files
vmap_file=open('versionmap.dat','w')
vmap_file.write('vnum\tL\tf\tphi\n')

task_file=open('run_210807_fPhi_PhaseDiag.task','w')
template_file='run_210807_fPhi_PhaseDiag.template'
template_contents=open(template_file,'r').read()


L = [4,5,6,7,8,9,10]
J=1.0
F = [0.1,0.2,0.4,0.6,0.8, 1.0, 1.1, 0.3]
phi_list = numpy.arange(0,4.1,0.1)
Nt=41
T = 1
Rmaxlist= [1,1,1,  1,1,4,6]
vnum=0
vnumt=0

for f in F:
        for phi in phi_list:
                phi = numpy.round(phi, 4)
                for l, Rmax in map(None, L, Rmaxlist):
                        for R in range(Rmax):
                                qsub_file=template_file.replace('.template','_'+str(vnum)+'.qsub')
                                fout=open(qsub_file,'w')
                                contents=template_contents.replace('###',str(vnum))
                                contents=contents.replace('*project*',project_name)
                                # *LLL*   *JJJ*   *fff*  *phi*  *Nt*  *TTT*
                                contents=contents.replace('*LLL*',str(l))
                                contents=contents.replace('*JJJ*',str(J))
                                contents=contents.replace('*fff*',str(f))
                                contents=contents.replace('*phi*',str(phi))
                                contents=contents.replace('*Nt*',str(Nt))
                                contents=contents.replace('*TTT*',str(T))
                                vmap_file.write(str(vnum)+'\t'+str(l)+'\t'+str(f)+'\t'+str(phi)+'\t'+str(Nt)+'\t'+str(T)+'\n')
                                if  l==9 and vnum%4==0:
                                        #if phi in [2.0,3.0,4.0]:
                                        #if f in [0.2, 0.4, 1.1] and 2.1<phi<4.0:
                                        #if phi in [4.0]:
                                        if f == 0.1 and 2.0<phi<4.0:
                                                print 'vnum: ', vnum,', J:', J,', L=', l,', f: ',f,', phi:', phi, ':  ',
                                                print 'bash run_210807_fPhi_PhaseDiag_'+str(vnum)+'.qsub'
                                                #print numpy.arange(-0.45*l,(d-1)*Omega*0.5,((d-1)*Omega*0.5-(-0.45*l))/12)
                                                task_file.write('bash run_210807_fPhi_PhaseDiag_'+str(vnum)+'.qsub\n')
                                                vnumt+=1
                                fout.write(contents)
                                fout.close()
                                vnum+=1




n_cpu_prtsk = 5
n_cores=(vnumt+0)
n_nodes=int(numpy.ceil(float(n_cpu_prtsk*vnumt)/partition_info[1]))
print n_cpu_prtsk, vnumt, n_nodes


for j in range(vnumt,n_cores):
        task_file.write('echo "Fake run"\n')

# Finally output sbatch file
contents=open('run_210807_fPhi_PhaseDiag.sbatch.template','r').read()
contents=contents.replace('*JobDetails*',str(' l==9 and f==0.1 and and 2.0<phi<4.0'))
contents=contents.replace('*nnn*',str(n_cores)) # The total number of processors
contents=contents.replace('*NNN*',str(n_nodes)) # The total number of nodes
contents=contents.replace('*ttt*',time_str) # The wall clock time per processor
contents=contents.replace('*ccc*',str(n_cpu_prtsk)) # Number of cpus per task
contents=contents.replace('*partition*',partition_info[0]) # Partition
contents=contents.replace('*myemail*',myemail) # Email
contents=contents.replace('*charge*',charge_str) # Project to charge to
open('run_210807_fPhi_PhaseDiag.sbatch','w').write(contents)





