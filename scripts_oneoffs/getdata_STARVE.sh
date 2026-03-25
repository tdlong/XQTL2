#!/bin/bash
#SBATCH --job-name=Demux
#SBATCH -A tdlong_lab
#SBATCH -p standard          

wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P01-TCTCATGA-TCCACGTG-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P01-TCTCATGA-TCCACGTG-READ2-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P02-CGAGGCCA-CCGGTAAC-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P02-CGAGGCCA-CCGGTAAC-READ2-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P03-TTCACGAG-GTAACAAT-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P03-TTCACGAG-GTAACAAT-READ2-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P04-GCGTGGAT-CATTGGTC-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P04-GCGTGGAT-CATTGGTC-READ2-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P05-TCTGGTAT-GTAAGCAA-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P05-TCTGGTAT-GTAAGCAA-READ2-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P06-CATTAGTG-TATGTAGT-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P06-CATTAGTG-TATGTAGT-READ2-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P07-ACGGTCAG-AACGAGGC-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P07-ACGGTCAG-AACGAGGC-READ2-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P08-GGCAAGCC-CGGATGCT-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P08-GGCAAGCC-CGGATGCT-READ2-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P09-TGTCGCTG-AGTCAGAC-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P09-TGTCGCTG-AGTCAGAC-READ2-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P10-ACCGTTAC-TCGCTATG-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P10-ACCGTTAC-TCGCTATG-READ2-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P11-TATGCCTT-TAATGTGT-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P11-TATGCCTT-TAATGTGT-READ2-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P12-ACTGGATC-CTAGCGGC-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P12-ACTGGATC-CTAGCGGC-READ2-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P13-TGGTACCT-AGTACTCA-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P13-TGGTACCT-AGTACTCA-READ2-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P14-TTGGAATT-GTATTGAC-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P14-TTGGAATT-GTATTGAC-READ2-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P15-CCTCTACA-AGGAGGTA-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P15-CCTCTACA-AGGAGGTA-READ2-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P16-GGAGCGTG-ACTTACGG-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P16-GGAGCGTG-ACTTACGG-READ2-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P17-GTCCGTAA-AAGATACA-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P17-GTCCGTAA-AAGATACA-READ2-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P19-TCAGAAGG-TATGATGG-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P19-TCAGAAGG-TATGATGG-READ2-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P20-GCGTTGGT-GGAAGTAT-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P20-GCGTTGGT-GGAAGTAT-READ2-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P21-ACATATCC-ATTGCACA-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P21-ACATATCC-ATTGCACA-READ2-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P22-TCATAGAT-CACCTTAA-READ1-Sequences.txt.gz
wget https://hts.igb.uci.edu/tdlong25032570/xR036-L2-G1-P22-TCATAGAT-CACCTTAA-READ2-Sequences.txt.gz
