#PBS -l select=8:ncpus=44
#PBS -l walltime=01:00:00
#PBS -j oe
#PBS -o D-reference.out

# Run the default version for the plane communication code
cd $PBS_O_WORKDIR
NPROCS=32

for NL in 1 2 4 8; do
    echo ---------
    echo Running on ${NL} nodes with ${NPROCS} cores per node
    aprun -n $((${NL} * ${NPROCS})) -N ${NPROCS} ./ft.D.x
    echo ---------
done
