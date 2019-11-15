#!/usr/bin/env bash

node_counts=(1 2 4 8 16 32 64 128 256 512)
files=(chapel.ft_transposed.D ref.mpi.ft.D ref.upc.ft-2d-upc.fftw3.DD)
files=(chapel.ft_transposed.E ref.mpi.ft.E ref.upc.ft-2d-upc.fftw3.DD8)
files=(chapel.ft_transposed.F ref.mpi.ft.F ref.upc.ft-2d-upc.fftw3.DD64)
for f in "${files[@]}"; do 
  for i in "${node_counts[@]}"; do
    out_file="$f.$i.out"
    if  [ -f "$out_file" ]; then 
      echo "$out_file"
      cat $out_file 2>/dev/null | grep "Elapsed time :"
      cat $out_file 2>/dev/null | grep "Time in seconds ="
      cat $out_file 2>/dev/null | grep "Total running time is"
      cat $out_file 2>/dev/null | grep "MFLOPS"
      cat $out_file 2>/dev/null | grep "Mop/s total"
      cat $out_file 2>/dev/null | grep "Mflops"
    echo ""
    fi
  done
  echo ""
  echo ""
done
