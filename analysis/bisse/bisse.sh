# load modules
module load boost/1.80.0-rnozcm7 
module load cmake/3.25.1-bc3eqck
module load autoconf/2.69-aorr3ao
module load automake/1.16.5-rpnnj5a

# add rb to the path
export PATH=$PATH:/work/LAS/phylo-lab/petrucci/revbayes/projects/cmake

# go to the correct directory
cd /work/LAS/phylo-lab/petrucci/ssetests_chap1/

# get the rep value
rep=$SLURM_ARRAY_TASK_ID

# create a file to hold the rep
touch aux/aux_${1}_$rep.Rev

# echo the definitions on it
printf "rep <- " >> aux/aux_${1}_$rep.Rev
echo $rep >> aux/aux_${1}_$rep.Rev

# source it, the parameter combination, and the actual script
timeout 24h rb aux/aux_${1}_$rep.Rev analysis/bisse/refs/refs_${1}.Rev analysis/bisse/master.Rev

# timeout and requeue if it takes more than 24h, MCMC is probably stuck
if [[ $? == 124 ]]; then 
  scontrol requeue ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
fi

# remove file
rm aux/aux_${1}_$rep.Rev
