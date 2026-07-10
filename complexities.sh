#!/bin/bash
ulimit -n 2048

for depth in 5 10 20 30; do
  for cplx in 0.1 0.2 0.4 0.6 0.8 1.0; do
    echo "Running depth=$depth cplx=$cplx"
    bash 05delPipeCoalComplex.sh \
      -g 1000000 \
      -t 0.01 \
      --depth $depth \
      --cplx $cplx \
      --pref G1000000P2T0.01D${depth}C${cplx} \
      -d ~/Desktop/The_University_of_Edinburgh/Sem2/Dissertation/Tetmer/CSKS_simulations/tmp
  done
done
