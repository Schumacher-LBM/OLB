#!/usr/bin/env bash
#2023 Nicolas Hafen

#Load .shellSupport for e.g. print command
. ~/.shellSupport/bash_main

#Run settings
runName='bifurcation3d'
pathInfo="tmp/batchRun/"
namePatch="patch.diff"
nameRunDoc="runInfo"
numOfThreads=8                  #Standalone (no slurm) only
executableName="bifurcation3d"  #Standalone (no slurm) only

#Parameters
resolutions=( 19 )
numberOfParticles=( 1000 )

## ==================== RUN ====================

#environment variables: https://slurm.schedmd.com/sbatch.html

#Check if standalone
if [ -z ${EXECUTABLE+x} ]; then 
  export SLURM_CLUSTER_NAME=$(hostname)
  export EXECUTABLE=$executableName
  export SLURM_NTASKS=$numOfThreads
  export SLURM_JOB_CPUS_PER_NODE=0 #TODO: fix
  export OMP_NUM_THREADS=$numOfThreads
  export PREFIX="-np $numOfThreads"
  export SLURM_SUBMIT_DIR=$(pwd)
  export MODULE_mpi=$(mpirun --version | head -n 1)
  export MODULE_compiler=$(gcc --version | head -n 1)
  #Create tmp, if not existing
  if [ ! -d "tmp" ]; then
    mkdir tmp
  fi
fi

#Create tmp/info, if not existing
if [ ! -d $pathInfo ]; then
  mkdir -p $pathInfo
fi

#Write patch (difference to current commit)
git diff > $pathInfo$namePatch

#Create runInfo
pathDoc=$pathInfo$nameRunDoc
startScriptTime=$(date +"%d.%m.%Y - %T")
if $rebuild ; then
  git log | head -n 3 > $pathDoc
  print "Starting run:     "$runName"\n" $pathDoc
else
  print "Continuing run:   "$runName"\n" $pathDoc
fi
print "Executable:       ${EXECUTABLE}\n" $pathDoc
print "Prefix:           ${PREFIX}\n" $pathDoc
print "Tasks:            ${SLURM_NTASKS}\n" $pathDoc
print "CPUs avail node:  ${SLURM_JOB_CPUS_PER_NODE}\n" $pathDoc
print "Threads:          ${OMP_NUM_THREADS}\n" $pathDoc
print "Directory:        ${SLURM_SUBMIT_DIR}\n" $pathDoc
print "Module compiler:  ${MODULE_compiler}\n" $pathDoc
print "Module mpi:       ${MODULE_mpi}\n" $pathDoc
print "Hostname:         $(hostname)\n" $pathDoc
print "Clustername:      ${SLURM_CLUSTER_NAME}\n" $pathDoc
print "Start time: $startScriptTime\n" $pathDoc
print "Parameters:" $pathDoc
#List resolutions
print "\n resolutions:" $pathDoc
for paraA in "${resolutions[@]}"
do
  print " $paraA," $pathDoc
done
#List number of particles
print "\n numberOfParticles:" $pathDoc
for paraB in "${numberOfParticles[@]}"
do
  print " $paraB," $pathDoc
done

print "\n\n///////////// Run Info //////////////\n\n\n" $pathDoc

#Do actual run
for N in "${resolutions[@]}"
do
  #Clean field
  make clean
  make
  for noP in "${numberOfParticles[@]}"
  do
    startTime=$(date +"%d.%m.%Y - %T")
    print "\n[======================================================]\n" $pathDoc
    print "[Parameters: N="$N" noP="$noP"\n" $pathDoc
    print "[Start time: $startTime\n" $pathDoc
    print "[======================================================]\n\n" $pathDoc
    #### RUN #####

    mpirun $PREFIX $EXECUTABLE ${N} ${noP}

    #### RUN ####

    #Save specified data (copy, to provide access to bifurcation.pvd)
    cp -r tmp/vtkData tmp/vtkDataN${N}P${noP}

    #Print converter and end time
    cat tmp/converter.dat >> $pathDoc 
    endTime=$(date +"%d.%m.%Y - %T")
    print "Start time: $startTime\n" $pathDoc
    print "End time: $endTime\n" $pathDoc
    print "\n\n\n" $pathDoc
  done
done

endScriptTime=$(date +"%d.%m.%Y - %T")
print "\n\nScript summary:\n" $pathDoc
print "Start time: $startScriptTime\n" $pathDoc
print "End time: $endScriptTime\n" $pathDoc

exit 0






