#!/bin/bash                                                                                             
#SBATCH -p sched_mit_g4nier                                                                             
#SBATCH -t 5-00:00:00    
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20                                                                                 
#SBATCH -J RAX4G_Sjob   
#SBATCH -o RAX4G_Sjob.out                                                                                         
#SBATCH --array=1-6


##add the modules
. /etc/profile.d/modules.sh
module add engaging/RAxML/8.2.9

##gets my array id, which is needed for use below. this will be, i think, a number like 1,2,3 etc
MY_ARRAY_ID=$SLURM_ARRAY_TASK_ID
echo $MY_ARRAY_ID                                                                 
THE_INDEX=~/4G_S/4G_S_Corr.txt
THE_INPUT_FILE=$( cat $THE_INDEX | grep "^$MY_ARRAY_ID " | awk '{print $2}' )
echo $THE_INPUT_FILE
RAXOUT=${THE_INPUT_FILE%%.*}

raxmlHPC-PTHREADS-AVX -T 20 -f a -m PROTGAMMALGF -p 12345 -x 12345 -#100 -n $RAXOUT -s $THE_INPUT_FILE

exit