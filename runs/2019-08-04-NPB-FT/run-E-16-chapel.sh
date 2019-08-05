#PBS -l select=16:ncpus=44
#PBS -l walltime=00:15:00
#PBS -j oe
#PBS -o E-16-chapel.out

# Run the default version for the plane communication code
cd $PBS_O_WORKDIR
NPROCS=32
NL=16
./ft -nl ${NL} --NPBClass=NPB.E
