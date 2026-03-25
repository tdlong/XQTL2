#!/bin/bash
#SBATCH --job-name=Demux
#SBATCH -A tdlong_lab
#SBATCH -p standard          

wget https://hts.igb.uci.edu/tdlong25042147/xR039-L8-G1-P4-CCACCAGGCA-ATTCCATAAG-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25042147/xR039-L8-G1-P4-CCACCAGGCA-ATTCCATAAG-READ2-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25042147/xR039-L8-G1-P5-GGTTGCGAGG-TTGCTCTATT-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25042147/xR039-L8-G1-P5-GGTTGCGAGG-TTGCTCTATT-READ2-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25042147/xR039-L8-G1-P6-ACAGTGTATG-CCGTATGTTC-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25042147/xR039-L8-G1-P6-ACAGTGTATG-CCGTATGTTC-READ2-Sequences.txt.gz

