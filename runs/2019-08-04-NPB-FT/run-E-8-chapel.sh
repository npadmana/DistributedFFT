#PBS -l select=8:ncpus=44
#PBS -l walltime=01:00:00
#PBS -j oe
#PBS -o E-8-chapel.out

# Run the default version for the plane communication code
cd $PBS_O_WORKDIR
NPROCS=32
NL=8
./ft -nl ${NL} --NPBClass=NPB.E
