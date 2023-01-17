# With this command line we simulate 3 times the same data set (since we use the same seed), but we have different output formats:
../../PolyMoSim-v1.1.3 -s 127 -m Modelfile_basic.txt -t tree1.txt -l log.log -n 1 -f site_pattern_freq_absolute  1> site_pattern_freq_abs.txt
../../PolyMoSim-v1.1.3 -s 127 -m Modelfile_basic.txt -t tree1.txt -l log.log -n 1 -f site_pattern_freq_relative  1> site_pattern_freq_rel.txt
../../PolyMoSim-v1.1.3 -s 127 -m Modelfile_basic.txt -t tree1.txt -l log.log -n 1 -f fasta                       1> sim-sequences.fas

../../PolyMoSim-v1.1.3 -s 127 -m Modelfile_basic.txt -t tree1.txt -l log.log -n 1 -f site_pattern_freq_absolute_fill  1> site_pattern_freq_abs_fill.txt
../../PolyMoSim-v1.1.3 -s 127 -m Modelfile_basic.txt -t tree1.txt -l log.log -n 1 -f site_pattern_freq_relative_fill  1> site_pattern_freq_rel_fill.txt

