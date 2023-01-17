../../PolyMoSim-v1.1.4 -s 11 -m Modelfile_basic.txt -t tree2.txt -l log.log --print_ancestral_seq anc.fas -n 1 --print_siterate_histogram hist.txt --print_siterate_data sitedata.txt -f fasta 1>sim-sequences.fas 2> sim-normal.log

../../PolyMoSim-v1.1.4 -s 11 -m Modelfile_basic.txt -t tree2.txt -l log.log --print_ancestral_seq anc.fas -n 1 --print_siterate_histogram hist.txt --print_siterate_data sitedata.txt -f site_pattern_freq_relative_fill 1> rel_fill_pattern_freq.txt 2> sim-rel-fill-freq.log

../../PolyMoSim-v1.1.4 -s 11 -m Modelfile_basic.txt -t tree2.txt -l log.log --print_ancestral_seq anc.fas -n 1 --print_siterate_histogram hist.txt --print_siterate_data sitedata.txt -f site_pattern_freq_relative 1> rel_pattern_freq.txt 2> sim-rel-freq.log