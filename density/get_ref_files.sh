pdbfixer --pdbid=9cbu --output=9cbu.pdb
sed -i 's/OP1/O1P/g' 9cbu.pdb
sed -i 's/OP2/O2P/g' 9cbu.pdb
pdbfixer --pdbid=7ez0 --output=7ez0.pdb
curl -LO https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-42499/map/emd_42499.map.gz
gzip -d emd_42499.map.gz
curl -LO https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-31385/map/emd_31385.map.gz
gzip -d emd_31385.map.gz

