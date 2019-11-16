echo "In run.sh script!"

function run_ref_mpi_perf {
  base=$(basename $2)
  out_file="ref.mpi.$base.$3.$1.out"
  echo "Running MPI $base size $3 with $1 nodes!"
  NPROCS=32
  aprun -n $((${1} * ${NPROCS})) -N ${NPROCS} ./$2.$3.x >> $out_file 2>&1
  echo "" >> $out_file

  echo -n "    " && date
}

function run_ref_upc_perf {
  export XT_SYMMETRIC_HEAP_SIZE=1024M
  base=$(basename $2)
  out_file="ref.upc.$base.$3.$1.out"
  echo "Running UPC $base size $3 with $1 nodes!"
  NPROCS=32
  aprun -n $((${1} * ${NPROCS})) -N ${NPROCS} ./$2.$3 ${NPROCS} ${1} >> $out_file 2>&1
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
node_counts=(1 2 4 8 16 32 64 128 256 512)
for i in {1..3}; do
  for nodes in "${node_counts[@]}"; do 
    run_perf     $nodes DistributedFFT/target/example/NPB-FT/ft_transposed D
    run_ref_mpi_perf $nodes NPB3.4-MPI/bin/ft D
    run_ref_upc_perf $nodes upc/ft-2d-upc.fftw3 DD
  done
done

# size E
node_counts=(8 16 32 64 128 256 512)
for i in {1..2}; do
  for nodes in "${node_counts[@]}"; do 
    run_perf         $nodes DistributedFFT/target/example/NPB-FT/ft_transposed E
    run_ref_mpi_perf $nodes NPB3.4-MPI/bin/ft E
    run_ref_upc_perf $nodes upc/ft-2d-upc.fftw3 DD8
  done
done

# size F
node_counts=(64 128 256 512)
for i in {1..1}; do
  for nodes in "${node_counts[@]}"; do 
    run_perf         $nodes DistributedFFT/target/example/NPB-FT/ft_transposed F
    run_ref_mpi_perf $nodes NPB3.4-MPI/bin/ft F
    run_ref_upc_perf $nodes upc/ft-2d-upc.fftw3 DD64
  done
done

date

cd $workdir

echo "DONE"
