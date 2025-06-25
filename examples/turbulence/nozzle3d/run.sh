#!/usr/bin/env sh
#
# This is a basic example script on how to utilize heteogeneous
# processing in OpenLB on a shared memory system. The underlying
# domain decompositions are generated using the
# `script/balancing/examples/balance_nozzle.py` notebook.
#
# Case must be built with CPU_SIMD and GPU_CUDA in hybrid mode,
# mixed config recommended. See `env-heterogeneity` in `flake.nix`.
#
# This is experimental. Only use if you have HPC experience.
# When in doubt ignore this file and use homogeneous MPI-only
# execution.

export N=9
export OLB="./nozzle3d --resolution ${N} --load-decomposition-from decomposition/n${N}.xml --heterogeneous --max-phys-t 100"

# adapt to number of NVIDIA GPUs in system
export N_GPU=2
# adapt to number of CPU cores in system
export N_CPU_CORE=16

mpirun --bind-to none \
	  -np $N_GPU -x OMP_NUM_THREADS=1           bash -c 'export CUDA_VISIBLE_DEVICES=${OMPI_COMM_WORLD_LOCAL_RANK}; ${OLB}' \
	: -np 1      -x OMP_NUM_THREADS=$N_CPU_CORE bash -c 'export CUDA_VISIBLE_DEVICES=;                              ${OLB}'
