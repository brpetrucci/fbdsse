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
timeout 121h rb aux/aux_${1}_$rep.Rev analysis/hisse/refs/refs_${1}.Rev analysis/hisse/master.Rev

# remove file
rm aux/aux_${1}_$rep.Rev
