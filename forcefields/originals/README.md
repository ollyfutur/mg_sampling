# Sources

- amber14sb.ff:
    - Amber14SB force field with ff99bsc0_chiOL3 for RNA, ported to Gromacs
    - [https://ftp.gromacs.org/contrib/forcefields/amber14sb.ff.tar.gz]
- Magnesium-FFs: 
    - microMg and nanoMg parameters for TIP3P, parmBSC0chiOL3 RNA, [Mamatkulov-Schwierz](https://doi.org/10.1063/1.5017694) ions
    - [https://github.com/bio-phys/Magnesium-FFs/archive/refs/heads/main.zip]
- MgFF_OPC: https://github.com/bio-phys/MgFF_OPC/archive/refs/heads/main.zip
    - microMg and nanoMg parameters for OPC (developped on SPC/E and shown transferable)
- Force-fields-for-metal-cations: https://github.com/bio-phys/Force-fields-for-metal-cations/archive/refs/heads/master.zip
    - parameters for the metal cations Li+, Na+, K+, Cs+, Mg2+, Ca2+, Sr2+, and Ba2+ in combination with the TIP3P
- amber14sb_OL15_O8scaled_v2.ff
    - amberOL15-ECC (Electronic continuum correction) forcefield provided by Duboué-Dijon
    - to be used with rescaled ionic charges (×0.8) (ions.ff) and SPC/E
- ions.ff
    - parameters for ions provided by Duboué-Dijon for use with amberOL15-ECC
- OPC:
    - OPC water model
    - adapted from [https://github.com/intbio/gromacs_ff/blob/2543183b3e68c9f504f5247e092e15317714b4a5/amber14sb_parmbsc1_opc_lmi.ff/opc.itp]
