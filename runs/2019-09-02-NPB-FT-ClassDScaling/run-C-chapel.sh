#PBS -l select=32:ncpus=44
#PBS -l walltime=01:00:00
#PBS -j oe
#PBS -o C-chapel.out

# Run the default version for the plane communication code
cd $PBS_O_WORKDIR
export TARGET=../../target/example/NPB-FT

chpl --version
module list

echo Running cases now....

for NL in 4 8 16 32; do
    echo ---------
    echo Running on ${NL} nodes
    ${TARGET}/ft_transposed -nl ${NL} --NPBClass=NPB.C
    echo ---------
done

