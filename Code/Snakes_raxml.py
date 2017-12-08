#!/usr/bin/python

# created abigailc@Artemis on november 8 2016

# this script will take a folder, and run all it's contents though muscle and raxml on the engaging cluster.

#set these yourself
ssh_inst = "ssh -l abigailc -i /mnt/c/Users/acaro/GradSchool/TreeOfLife/Snakes/id_rsa eofe5.mit.edu"
clus_head = "abigailc@eofe5.mit.edu:/home/abigailc/"

#ssh_inst = "ssh dgruen eofe5.mit.edu"
#clus_head = "dgruen@eofe5.mit.edu:/home/dgruen/"

#imports
import sys
import argparse
import os
import re
import time

#makes dir if need be
def check_directory_existance(prefix, ssh_inst):
    import os
    print("checking dirs")
    os.system(ssh_inst+" \'mkdir "+prefix+"\'")


def gen_correlate_file(list_of_input_files, corr_file):
    #this should be in form
    #1 name
    #2 name
    #3 name
    #requires 1: list of files 2. name for corr_file.
    i = 0
    with open(corr_file, "w") as corr:
        for item in list_of_input_files:
            i += 1
            #make sure the \n spacing works correctly.
            corr.write(str(i)+" "+item+"\n")
    return corr_file



def gen_raxml_array_script(scriptfile, indexname, n, Jobname, AA_MODEL):
    #currently assuming you are running in the dir that files are in and should be returned to.
    #thats it. just print the script. return its filename, which will need to be added to list of things to be moved to the cluster.
##example script
    a =  """#!/bin/bash                                                                                             
#SBATCH -p sched_mit_g4nier                                                                             
#SBATCH -t 5-00:00:00    
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20                                                                                 
#SBATCH -J RAX"""+Jobname+"""   
#SBATCH -o RAX"""+Jobname+""".out                                                                                         
#SBATCH --array=1-"""+n+"""


##add the modules
. /etc/profile.d/modules.sh
module add engaging/RAxML/8.2.9

##gets my array id, which is needed for use below. this will be, i think, a number like 1,2,3 etc
MY_ARRAY_ID=$SLURM_ARRAY_TASK_ID
echo $MY_ARRAY_ID                                                                 
THE_INDEX="""+indexname+"""
THE_INPUT_FILE=$( cat $THE_INDEX | grep "^$MY_ARRAY_ID " | awk '{print $2}' )
echo $THE_INPUT_FILE
RAXOUT=${THE_INPUT_FILE%%.*}

raxmlHPC-PTHREADS-AVX -T 20 -f a -m """+AA_MODEL+""" -p 12345 -x 12345 -#100 -n $RAXOUT -s $THE_INPUT_FILE

exit"""
    with open(scriptfile, "w") as script:
        script.write(a)
    return scriptfile

    

#does everything
def raxml_array_on_cluster(input_list, prefix, AA_MODEL):
    #this creates dir you will use on the cluster.
    finished_list = []
    remove_list = []
    for item in input_list:
        exts = item.split(".")
        tail = exts[0]
        head = "RAxML_bipartitions."
        new = head+tail
        finished_list.append(new)
    check_directory_existance(prefix, ssh_inst)
    clus_path = "/"+prefix
    a = gen_raxml_array_script(prefix+"_Sc.sh", "~"+clus_path+"/"+prefix+"_Corr.txt", str(len(input_list)), prefix+"job", AA_MODEL)
    b = gen_correlate_file(input_list, prefix+"_Corr.txt")
    input_list.append(a)
    input_list.append(b)
    direct = os.getcwd()
    move_to_cluster(input_list, clus_path)
    print("everything should be generated and on the cluster")
    n = str(len(input_list))
    print("there are "+n+" files that should be RAXML-ing right now")
    os.system(ssh_inst+" 'cd ~/"+prefix+";echo $PWD;sbatch "+a+"'")
    finished = "start"
    #to see if the run is complete, see if each new file has been generated. check every 5 minutes for muscle.
    #initialize list of things to move home. initially equal to aligned_list.
    movehome = []
    time.sleep(180)
   
    for i in finished_list:
        movehome.append(i)
    while finished is not True:
        #try and move each home.
        for filename in movehome:
            os.system("scp "+clus_head[:-1]+clus_path+"/"+filename+" "+direct)
        finished = "yes"
        for item in finished_list:
            #see if it got moved home.
            exists = os.path.isfile(item)
            if exists is True:
                if item in movehome:
                    movehome.remove(item)
            else:
                finished = False
                print(item+" not found")
        if finished == "yes":
            print("Should be done!")
            finished = True
        else:
            #wait five minutes and then try again.
            print("checking.... some alignment outputs do not exist yet. sleeping for 10 minutes.")
            time.sleep(600)
    #files should        
    print("Your files have been aligned! They are located at "+direct)
    return finished_list



#moves your .fasta and script/corr to cluster
def move_to_cluster(list_of_files, clus_path):
    #requires os.system
    #requires scp(?)
    #do the thing
    for item in list_of_files:
        os.system("scp "+item+" "+clus_head+clus_path)
    print("Finished moving files to cluster in place:"+clus_path)

# parser
#removed -- check Sinker.py for changes.



