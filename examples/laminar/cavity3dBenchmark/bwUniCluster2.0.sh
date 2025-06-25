#!/usr/bin/env bash
#SBATCH --partition="dev_gpu_4"  ####--> REF:https://wiki.bwhpc.de/e/BwUniCluster2.0/Batch_Queues
#SBATCH --gres=gpu:4
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --time=00:30:00
#SBATCH --job-name=cavity3dBenchmark_GPU
#SBATCH --error=job.\%J.err
#SBATCH --output=job.\%J.out

echo "--------------- module purge and load ---------------"
module purge 
module load toolkit/nvidia-hpc-sdk/23.9
echo "--------------- end module purge and load ---------------"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.
#export WORKSPACE=/pfs/work7/workspace/scratch/ve6116-CFDSLab_P2_test

echo "--------------- Builing and nvidia-smi ---------------"
make clean-core;make core;make clean;make;nvidia-smi
echo "--------------- End builing  ---------------"

echo "--------------- Symbolic link   ---------------"
if [ -d tmp ]; then
   if [ -L tmp ]; then
      unlink tmp
   else
      rm -rf tmp
   fi
fi
mkdir $WORKSPACE/job_${SLURM_JOB_ID}
ln -s $WORKSPACE/job_${SLURM_JOB_ID} tmp
echo "--------------- End symbolic link   ---------------"

for n in 1 2 3 4 ; do
    echo "--------------- ${n} GPU are running---------------"
    mpirun  --mca btl_openib_warn_default_gid_prefix 0  -np ${n} bash -c 'export CUDA_VISIBLE_DEVICES=${OMPI_COMM_WORLD_LOCAL_RANK}; ./cavity3d --size 200 --no-results'  2>&1 |tee outlog${n}.dat
    ####  -mca pml ucx --mca btl_openib_warn_default_gid_prefix 0 is in https://wiki.bwhpc.de/e/BwUniCluster2.0/Slurm
    echo "--------------- ${n} GPU end---------------"
done


###see in https://wiki.bwhpc.de/e/BwUniCluster2.0/Slurm#Slurm_Commands_.28excerpt.29 how to submit a job
