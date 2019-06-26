#PBS -l select=16:ncpus=44
#PBS -l walltime=01:00:00
#PBS -j oe
#PBS -o plane-communication

# Run the default version for the plane communication code
cd $PBS_O_WORKDIR

for NG in 128 512 1024 1536 2048 2560; do
    ../../target/example/plane -nl 16 --Ng=$NG
done
