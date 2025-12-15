#!/bin/bash

for id in 9CBU 9CBW 9CBX 9CBY 7EZ0; do
   curl -LO https://files.rcsb.org/download/$id.pdb
   pdbfixer $id.pdb --add-atoms=heavy --keep-heterogens=all --output=${id}_fixed.pdb
done

