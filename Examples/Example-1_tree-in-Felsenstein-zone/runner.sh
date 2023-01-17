# The following command line will run the simulation of this example:
../../PolyMoSim-v1.1.4 -s 127 -m Modelfile_basic.txt -t tree1.txt -l log.log --print_ancestral_seq anc.fas -n 1 --print_siterate_histogram hist.txt --print_siterate_data sitedata.txt  1>sim-sequences_1k.nex

../../PolyMoSim-v1.1.4 -s 127 -m Modelfile_basic.txt -t tree1.txt -l log.log -n 1 -f site_pattern_freq_absolute  1> site_pattern_freq_abs_1k.txt
../../PolyMoSim-v1.1.4 -s 127 -m Modelfile_basic.txt -t tree1.txt -l log.log -n 1 -f site_pattern_freq_relative  1> site_pattern_freq_rel_1k.txt
../../PolyMoSim-v1.1.4 -s 127 -m Modelfile_basic.txt -t tree1.txt -l log.log -n 1 -f site_pattern_freq_relative_fill  1> site_pattern_freq_rel_fill_1k.txt
../../PolyMoSim-v1.1.4 -s 127 -m Modelfile_basic.txt -t tree1.txt -l log.log -n 1 -f fasta                       1> sim-sequences_1k.fas

