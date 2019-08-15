#!/usr/bin/env bash

node_counts=(1 2 4 8 16 32 64 128 256)
files=(chapel.ft_transposed.D ref.ft.D)
files=(chapel.ft_transposed.E ref.ft.E)
#files=(chapel.ft_transposed.F ref.ft.F)
for f in "${files[@]}"; do 
  for i in "${node_counts[@]}"; do
    out_file="$f.$i.out"
    if  [ -f "$out_file" ]; then 
      echo "$out_file"
      cat $out_file 2>/dev/null | grep "Elapsed time :"
      cat $out_file 2>/dev/null | grep "Time in seconds ="
      #cat $out_file 2>/dev/null | grep "MFLOPS"
      #cat $out_file 2>/dev/null | grep "Mop/s total"
    echo ""
    fi
  done
done
