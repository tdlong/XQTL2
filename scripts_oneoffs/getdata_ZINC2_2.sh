#!/bin/bash
#SBATCH --job-name=Demux
#SBATCH -A tdlong_lab
#SBATCH -p standard          

wget https://hts.igb.uci.edu/tdlong25042147/xR039-L8-G1-P1-TAAGCATCCA-AATGGATTGA-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25042147/xR039-L8-G1-P1-TAAGCATCCA-AATGGATTGA-READ2-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25042147/xR039-L8-G1-P2-ACCACGACAT-CCGCATACGA-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25042147/xR039-L8-G1-P2-ACCACGACAT-CCGCATACGA-READ2-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25042147/xR039-L8-G1-P3-GCCGCACTCT-CGAGGTCGGA-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25042147/xR039-L8-G1-P3-GCCGCACTCT-CGAGGTCGGA-READ2-Sequences.txt.gz

