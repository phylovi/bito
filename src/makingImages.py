import bito
import os

# Set working directory to your libsbn directory location
os.chdir(os.path.join(os.environ["HOME"], 'matsen/libsbn'))

fasta_file = 'data/ds1-reduced-5.fasta'
newick_file = 'data/ds1-reduced-5.nwk'

inst = bito.gp_instance('_ignore/mmapped_plv.data')
inst.read_fasta_file(fasta_file)
inst.read_newick_file(newick_file)
inst.make_engine()
# parent1 = "00101"
# parent2 = "00010"
# child1 = "00100"
# child2 = "00001"
# inst.add_node_pair(parent1, parent2, child1, child2)
inst.subsplit_dag_to_dot('_ignore/out_ds1_5.dot', True)


inst2 = bito.gp_instance('_ignore/mmapped_plv.data')
inst2.read_fasta_file("data/7-taxon-slice-of-ds1.fasta")
inst2.read_newick_file("data/simplest-hybrid-marginal-all-trees.nwk")
inst2.make_engine()
# parent1 = "01100"
# parent2 = "00011"
# child1 = "01000"
# child2 = "00100"
# inst2.add_node_pair(parent1, parent2, child1, child2)
# inst2.subsplit_dag_to_dot('_ignore/out_prep_add.dot', True)
# parent1 = "10000"
# parent2 = "01100"
# child1 = "01000"
# child2 = "00100"
# inst2.add_node_pair(parent1, parent2, child1, child2)
inst2.subsplit_dag_to_dot('_ignore/out_ds1_7_all.dot', True)

inst3 = bito.gp_instance('_ignore/mmapped_plv.data')
inst3.read_fasta_file("data/7-taxon-slice-of-ds1.fasta")
inst3.read_newick_file("data/simplest-hybrid-marginal.nwk")
inst3.make_engine()
# parent1 = "01000"
# parent2 = "00111"
# child1 = "00101"
# child2 = "00010"
# inst3.add_node_pair(parent1, parent2, child1, child2)
inst3.subsplit_dag_to_dot('_ignore/out_ds1_7_simple.dot', True)

inst4 = bito.gp_instance('_ignore/mmapped_plv.data')
inst4.read_fasta_file("data/7-taxon-slice-of-ds1.fasta")
inst4.read_newick_file("data/second-simplest-hybrid-marginal.nwk")
inst4.make_engine()
# parent1 = "01000"
# parent2 = "00111"
# child1 = "00101"
# child2 = "00010"
# inst3.add_node_pair(parent1, parent2, child1, child2)
inst4.subsplit_dag_to_dot('_ignore/out_ds1_7_simple_2.dot', True)