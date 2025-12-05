#!/bin/bash

curl -LO https://files.rcsb.org/download/9CBU.pdb
curl -LO https://files.rcsb.org/download/9CBW.pdb
curl -LO https://files.rcsb.org/download/9CBX.pdb
curl -LO https://files.rcsb.org/download/9CBY.pdb
curl -LO https://files.rcsb.org/download/7EZ0.pdb
curl -LO https://files.wwpdb.org/pub/emdb/structures/EMD-42498/map/emd_42498.map.gz
curl -LO https://files.wwpdb.org/pub/emdb/structures/EMD-42498/other/emd_42498_half_map_1.map.gz
curl -LO https://files.wwpdb.org/pub/emdb/structures/EMD-42498/other/emd_42498_half_map_2.map.gz
curl -LO https://files.wwpdb.org/pub/emdb/structures/EMD-42499/map/emd_42499.map.gz
curl -LO https://files.wwpdb.org/pub/emdb/structures/EMD-42499/other/emd_42499_half_map_1.map.gz
curl -LO https://files.wwpdb.org/pub/emdb/structures/EMD-42499/other/emd_42499_half_map_2.map.gz
gzip -d emd_42498.map.gz emd_42498_half_map_1.map.gz emd_42498_half_map_2.map.gz emd_42499.map.gz emd_42499_half_map_1.map.gz emd_42499_half_map_2.map.gz
curl -LO https://files.wwpdb.org/pub/emdb/structures/EMD-31385/map/emd_31385.map.gz
gzip -d emd_31385.map.gz
for f in *.map; do
   mv $f ${f%.*}.mrc
done
