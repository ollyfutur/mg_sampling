#!/bin/bash

#dirnames="ob_mMg_opc ob_nMg_tip3p op_mMg_opc op_nMg_tip3p ob_amber_tip3p ob_mMg_tip3p op_amber_tip3p op_mMg_tip3p"
dirnames="op_mMg_tip3p"
for dirname in $dirnames; do
  (
  cd $dirname
  plumed sum_hills --hills HILLS
  )
done

