#!/bin/bash

#refs="MG NA"
refs="MG NA CL"
sels="RNA O_Water"
for ref in $refs; do
  for sel in $sels; do
      $GMX rdf -f 0_trajout.xtc -s ../rep00/run.tpr -n ../index.ndx -o 0_${ref}_${sel}_rdf.xvg -b 75000 -ref ${ref} -sel ${sel} &
  done
  wait
done
