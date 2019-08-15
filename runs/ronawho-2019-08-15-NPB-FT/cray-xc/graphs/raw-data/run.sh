echo "In run.sh script!"

function run_ref_perf {
  base=$(basename $2)
  out_file="ref.$base.$3.$1.out"
  echo "Running $2 size $3 with $1 nodes!"
  NPROCS=32
  aprun -n $((${1} * ${NPROCS})) -N ${NPROCS} ./$2.$3.x >> $out_file 2>&1
  echo "" >> $out_file

  echo -n "    " && date
}

function run_perf {
  base=$(basename $2)
  out_file="chapel.$base.$3.$1.out"
  echo "Running $base size $3 with $1 nodes!"
  ./$2 --NPBClass=NPB.$3 -nl $1 >> $out_file 2>&1
  echo "" >> $out_file

  echo -n "    " && date
}

date

chpl --version
module list 

# size D
node_counts=(1 2 4 8 16 32 64 128 256)
for i in {1..2}; do
  for nodes in "${node_counts[@]}"; do 
    run_perf     $nodes DistributedFFT/target/example/NPB-FT/ft_transposed D
    run_ref_perf $nodes NPB3.4-MPI/bin/ft D
  done
done

# size E
node_counts=(8 16 32 64 128 256)
for i in {1..2}; do
  for nodes in "${node_counts[@]}"; do 
    run_perf     $nodes DistributedFFT/target/example/NPB-FT/ft_transposed E
    run_ref_perf $nodes NPB3.4-MPI/bin/ft E
  done
done

# size F
node_counts=(64 128 256)
#node_counts=(512)
for i in {1..2}; do
  for nodes in "${node_counts[@]}"; do 
    run_perf     $nodes DistributedFFT/target/example/NPB-FT/ft_transposed F
    run_ref_perf $nodes NPB3.4-MPI/bin/ft F
  done
done

date

cd $workdir

echo "DONE"
