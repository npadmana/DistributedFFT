#PBS -l select=8:ncpus=44
#PBS -l walltime=01:00:00
#PBS -j oe
#PBS -o E-8-reference.out

# Run the default version for the plane communication code
cd $PBS_O_WORKDIR
NPROCS=32
NL=8
aprun -n $((${NL} * ${NPROCS})) -N ${NPROCS} ./ft.E.x
