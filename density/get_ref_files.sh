#!/bin/bash

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

# Con2-2.2A half maps of EMD-42499 (DasLab/Water-CryoEM-ribozyme), header set to volume (ispg=1)
con2=https://raw.githubusercontent.com/DasLab/Water-CryoEM-ribozyme/d470420b9dfa178876493c4e282611fecf86db1e/maps
curl -LO $con2/emd_42499_Con2-2.2A_half_A.mrc
curl -LO $con2/emd_42499_Con2-2.2A_half_B.mrc
python - <<EOF
import mrcfile
for f in ["emd_42499_Con2-2.2A_half_A.mrc", "emd_42499_Con2-2.2A_half_B.mrc"]:
    with mrcfile.open(f, mode="r+") as m:
        m.set_volume()
EOF
