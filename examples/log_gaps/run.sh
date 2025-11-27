#!/bin/bash

seed=$1

for(( method=1; method<=2; method++))
do
  for n in 100 150 200 250
  do
    R CMD BATCH "--args seed=$seed method=$method n=$n" log_gaps.R method$method\_n$n\_seed$seed.Rout &
  done
done

for n in 100 150 200 250
  do
    python3 log_gaps.py $seed 1 $n &>/dev/null &
done
