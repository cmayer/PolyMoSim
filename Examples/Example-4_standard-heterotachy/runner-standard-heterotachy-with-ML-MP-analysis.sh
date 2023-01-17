../../PolyMoSim-v1.1.2  -m Modelfile_basic.txt -t tree-4-standard-heterotachy.txt -l log.log --print_ancestral_seq anc.fas -n 1 -f phylip 1> sim-sequences.phy

## If you have raxml installed you can check the reconstruciton success:
#  In most cases both, the MP as well as the ML method will yield the wrong topology.

raxml -y -s sim-sequences.phy -m GTRGAMMA -n parsimony_tree -p 123
raxml -s sim-sequences.phy -m GTRGAMMA -n likelihood_tree -p 123
