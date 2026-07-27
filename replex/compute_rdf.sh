#!/bin/bash

#refs="MG NA"
refs="MG NA CL"
sels="RNA O_Water"
reps="0_trajout.xtc"
for ref in $refs; do
  for sel in $sels; do
    for rep in $reps; do
      $GMX rdf -f ${rep}/run.xtc -s ${rep}/run.tpr -n index.ndx -o ${rep}/${ref}_${sel}_rdf.xvg -b 50000 -ref ${ref} -sel ${sel} &
    done
    wait
  done
done
