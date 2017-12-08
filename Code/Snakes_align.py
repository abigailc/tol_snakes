#!/usr/bin/python

# created abigailc@Artemis on november 8 2016
# edited abigailc@Leviathan on december 3 2017 7pm

# ONLY WORKS on SLURM scheduling clusters
# this script will take a folder, and run all it's contents though muscle on the engaging cluster.

#set these yourself
ssh_inst = "ssh -l abigailc -i /mnt/c/Users/acaro/GradSchool/TreeOfLife/Snakes/id_rsa eofe5.mit.edu"
clus_head = "abigailc@eofe5.mit.edu:/home/abigailc/"

#imports
import sys
import argparse
import os
import re
import time


#functions

#makes dir if need be
def check_directory_existance(prefix, ssh_inst):
    import os
    print("checking dirs")
    os.system(ssh_inst+" \'mkdir "+prefix+"\'")

#makes a correlated file - number file
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

#makes a muscle script with correct name and correlates
def gen_muscle_array_script(scriptfile, indexname, n, Jobname):
    #currently assuming you are running in the dir that files are in and should be returned to.
    #thats it. just print the script. return its filename, which will need to be added to list of things to be moved to the cluster.
##example script
    a =  """#!/bin/bash                                                                                             
#SBATCH -p sched_mit_g4nier                                                                             
#SBATCH -t 5-00:00:00    
#SBATCH --nodes=1                                                                               
#SBATCH -J RAX"""+Jobname+"""   
#SBATCH -o RAX"""+Jobname+""".out                                                                                         
#SBATCH --array=1-"""+n+"""


##add the modules
. /etc/profile.d/modules.sh
module add engaging/muscle/3.8.31

##gets my array id, which is needed for use below. this will be, i think, a number like 1,2,3 etc
MY_ARRAY_ID=$SLURM_ARRAY_TASK_ID
echo $MY_ARRAY_ID                                                                 
THE_INDEX="""+indexname+"""
THE_INPUT_FILE=$( cat $THE_INDEX | grep "^$MY_ARRAY_ID " | awk '{print $2}' )
echo $THE_INPUT_FILE
TAIL=_Muscle.fasta
HEAD=${THE_INPUT_FILE%%.*}
MUSOUT=$HEAD$TAIL

muscle -in $THE_INPUT_FILE -out $MUSOUT

exit"""
    with open(scriptfile, "w") as script:
        script.write(a)
    return scriptfile

    

#does everything
def muscle_array_on_cluster(input_list, prefix):
    #this creates dir you will use on the cluster.
    aligned_list = []
    remove_list = []
    for item in input_list: #make list of aligned file names
        if "_Sc.sh" in item:
            continue
        if "_Corr.txt" in item:
            continue
        bits = item.split(".")
        head = bits[0]
        new = head+"_Muscle.fasta"
        aligned_list.append(new)
    check_directory_existance(prefix, ssh_inst) #on the cluster
    clus_path = "/"+prefix
    a = gen_muscle_array_script(prefix+"_Sc.sh", "~"+clus_path+"/"+prefix+"_Corr.txt", str(len(input_list)), prefix+"job")
    b = gen_correlate_file(input_list, prefix+"_Corr.txt")
    input_list.append(a)
    input_list.append(b)
    direct = os.getcwd()
    move_to_cluster(input_list, clus_path)
    print("everything should be generated and on the cluster")
    n = str(len(input_list))
    print("there are "+n+" files that should be aligning in muscle right now")
    os.system(ssh_inst+" 'cd ~/"+prefix+";echo $PWD;sbatch "+a+"'")
    finished = "start"
    #to see if the run is complete, see if each new file has been generated. check every 5 minutes for muscle.
    #initialize list of things to move home. initially equal to aligned_list.
    movehome = []
    time.sleep(180)
   
    for i in aligned_list:
        movehome.append(i)
    while finished is not True:
        #try and move each home.
        for filename in movehome:
            os.system("scp "+clus_head[:-1]+clus_path+"/"+filename+" "+direct)
        for item in aligned_list:
            #see if it got moved home.
            exists = os.path.isfile(item)
            if exists is True:
                if item in movehome:
                    movehome.remove(item)
                finished = "yes"
            else:
                finished = False
                print(item+" not found")
        if finished == "yes":
            print("Should be done!")
            finished = True
        else:
            #wait five minutes and then try again.
            print("checking.... some alignment outputs do not exist yet. sleeping for 5 minutes.")
            time.sleep(300)
    #files should        
    print("Your files have been aligned! They are located at "+direct)
    return aligned_list



#moves your .fasta and script/corr to cluster
def move_to_cluster(list_of_files, clus_path):
    #requires os.system
    #requires scp(?)
    #do the thing
    for item in list_of_files:
        os.system("scp "+item+" "+clus_head+clus_path)
    print("Finished moving files to cluster in place:"+clus_path)

# parser

if __name__ == "__main__":
    print("Running in terminal")  
    import sys
    import argparse
    import os
    import re
    import time
    parser = argparse.ArgumentParser(description="All")
    parser.add_argument("-p", "--projectname", action = "store", default = "Unnamed_Job", help="give a name for your script/job")
    parser.add_argument("-f", "--folder", action = "store", default = False, help="give a folder name. should contain ONLY .fasta files.")
    
    args = parser.parse_args()
    #create input file list

    Input_List = os.listdir(args.folder)
    #change the directory so we are within args.folder
    try:
        os.chdir(args.folder)
    except:
        print ("didn't change dir")
        raise SystemExit
    #run the thing
    muscle_array_on_cluster(Input_List, args.projectname)
    print("done")



