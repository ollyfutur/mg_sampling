#!/bin/bash

refs="MG NA"
sels="O_Water RNA"
reps="rep??"
for ref in $refs; do
  for sel in $sels; do
    for rep in $reps; do
      $GMX rdf -f ${rep}/run.xtc -s ${rep}/run.tpr -n index.ndx -o ${rep}/${ref}_${sel}_rdf.xvg -b 50000 -ref ${ref} -sel ${sel} &
    done
    wait
  done
done
