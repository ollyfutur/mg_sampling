### Generate index file
$GMX make_ndx -f run.tpr -o index.ndx << EOF
"SOL" & a OW
"SOL" & a HW1 HW2
a NA
a CL
q
EOF

sed -i "s/SOL_&_OW/O_Water/g" index.ndx
sed -i "s/SOL_&_HW1_HW2/H_Water/g" index.ndx

for ref in NA CL; do
  for sel in O_Water H_Water NA CL; do
    $GMX rdf -f run.xtc -s run.tpr -n index.ndx -dt 10 -ref ${ref} -sel ${sel} -o rdf_${ref}_${sel}.xvg
  done
done
