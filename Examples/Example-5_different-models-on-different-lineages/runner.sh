rm best-paup-MP-tree.tre
rm best-paup-ML-tree.tre
MultiMoSeqSim -m Modelfile5.txt -t tree5.txt -l log.log  -n 10 --post paup-block-ML_MP-analysis.txt   1>sim-sequences.nex

## If you have the good old paup installed you can run the following test:
## To run the test, remove the # symbols in the following lines.

#paup -n sim-sequences.nex
#grep "(C,D)" best-paup-MP-tree.tre | wc
#grep "(C,D)" best-paup-ML-tree.tre | wc

MultiMoSeqSim -m Modelfile5.txt -f phylip -t tree5.txt -l log2.log  -n 1  1>sim-sequences.phy
