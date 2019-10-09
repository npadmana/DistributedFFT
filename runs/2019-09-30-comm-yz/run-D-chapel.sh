#PBS -l select=32,place=scatter
#PBS -l walltime=01:00:00
#PBS -j oe

# Run the default version for the plane communication code
cd $PBS_O_WORKDIR
export TARGET=../../target/example/NPB-FT

chpl --version
module list

echo Running cases now....

for NL in 8 16 32; do
    echo ---------
    echo Running on ${NL} nodes
    ${TARGET}/ft_transposed -nl ${NL} --NPBClass=NPB.D --timeTrackFFT
    ${TARGET}/ft_transposed -nl ${NL} --NPBClass=NPB.D --timeTrackFFT --overSubscribe=false
    echo ---------
done

