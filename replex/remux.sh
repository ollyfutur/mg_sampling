#!/bin/bash

mkdir -p remuxed

(
  shopt -s extglob
  cd remuxed
  $GMX trjcat -demux ../rep??/bias_index.xvg -f ../rep??/run.xtc
)
