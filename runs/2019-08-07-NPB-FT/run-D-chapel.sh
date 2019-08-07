#PBS -l select=8:ncpus=44
#PBS -l walltime=01:00:00
#PBS -j oe
#PBS -o D-chapel.out

# Run the default version for the plane communication code
cd $PBS_O_WORKDIR
export TARGET=../../target/example/NPB-FT

##echo //////////////////////////////////
##echo Default affinity
##export QT_AFFINITY=
##echo QT_AFFINITY=${QT_AFFINITY}
##echo //////////////////////////////////
##for NL in 1 8; do
##    echo ---------
##    echo Running on ${NL} nodes 
##    ${TARGET}/ft -nl ${NL} --NPBClass=NPB.D
##    echo ---------
##done

echo //////////////////////////////////
echo No Affinity
export QT_AFFINITY=no
echo QT_AFFINITY=${QT_AFFINITY}
echo //////////////////////////////////
for NL in 1 8; do
    echo ---------
    echo Running on ${NL} nodes 
    ${TARGET}/ft -nl ${NL} --NPBClass=NPB.D
    echo ---------
done

##echo //////////////////////////////////
##echo OpenMP runtime, default affinity
##export QT_AFFINITY=
##echo QT_AFFINITY=${QT_AFFINITY}
##echo //////////////////////////////////
##export QT_AFFINITY=
##for NL in 1 8; do
##    echo ---------
##    echo Running on ${NL} nodes 
##    ${TARGET}/ft_omp -nl ${NL} --NPBClass=NPB.D
##    echo ---------
##done

echo //////////////////////////////////
echo OpenMP runtime, no affinity
export QT_AFFINITY=no
echo QT_AFFINITY=${QT_AFFINITY}
echo //////////////////////////////////
export QT_AFFINITY=no
for NL in 1 8; do
    echo ---------
    echo Running on ${NL} nodes 
    ${TARGET}/ft_omp -nl ${NL} --NPBClass=NPB.D
    echo ---------
done
