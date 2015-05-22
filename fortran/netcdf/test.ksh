#!/bin/ksh
#
##################################################
# For btc runs

# set the partition where the job will run
#SBATCH --partition=compute

# set the number of nodes
#SBATCH --nodes=1

# tasks per node
#SBATCH --ntasks-per-node=24

# set max wallclock time
#SBATCH --time=00:10:00

# set name of job
#SBATCH --job-name="pnetcdf"

# mail alert at start, end and abortion of execution
##SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=lars.buntemeyer@hzg.de

#
#SBATCH --account=ch0636

## output file
##SBATCH --output ./out.o
#
## error file
##SBATCH --error ./error.e
#
##################################################
#
ulimit -s unlimited

#rm testfile
#echo 'striping testfile'
#lfs setstripe -c -1 testfile
#lfs getstripe testfile

time srun nc4_pres_temp_4D_wr

exit

