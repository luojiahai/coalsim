import qiuyi_tree
import pprint
import os, shutil
import numpy as np


def main():
    shutil.rmtree('./output')
    os.mkdir('./output')
    os.mkdir('./output/newick_gene_subtrees')
    os.mkdir('./output/subtrees')

    qstree = qiuyi_tree.SpeciesTree(newick_path='data/tree_sample.txt')    # read newick species tree
    qstree.save_to_file(path='output/species_nodes_table.txt')

    qiuyi_tree.SpeciesTree.global_species_tree = qstree
    qiuyi_tree.SpeciesTree.lambda0 = np.random.gamma(shape=3, scale=0.1, size=len(qstree.leaves))

    print('\nsecies_tree ascii_art:')
    print(qstree.skbio_tree.ascii_art())
    print('\nspecies_nodes:')
    qstree.print_nodes()
    print('\ncoalescent:')
    coalescent_process, _ = qstree.coalescent(distance_above_root=10000)      # do coalescece based on the species tree
    print('\ncoalescent_process:')
    pprint.pprint(coalescent_process)
    print('\ntime_sequences:')
    time_sequences = qstree.time_sequences(coalescent_process=coalescent_process)       # convert to time sequence structure
    pprint.pprint(time_sequences)
    
    # print('\nTEST')
    # pprint.pprint(qstree.filter_coal_process(coalescent_process, '1*2*'))

    qgtree = qiuyi_tree.GeneTree(time_sequences=time_sequences, species_tree=qstree)        # construct newick coalescent tree

    qiuyi_tree.GeneTree.lambda_dup = np.random.gamma(shape=2, scale=0.1, size=len(qgtree.leaves))
    qiuyi_tree.GeneTree.lambda_loss = np.random.gamma(shape=2, scale=0.1, size=len(qgtree.leaves))
    qiuyi_tree.GeneTree.lambda_trans = np.random.gamma(shape=2, scale=0.1, size=len(qgtree.leaves))

    qgtree.save_to_file(path='output/gene_nodes_table.txt')
    f = open('output/newick_gene_subtrees/gene_tree.txt', 'w')
    f.write(str(qgtree.skbio_tree))
    f.close()
    print('\ngene_tree ascii_art:')
    print(qgtree.skbio_tree.ascii_art())
    print('\ngene_nodes:')
    qgtree.print_nodes()
    print('\ngene_tree dlt_process:')
    events = qgtree.dup_loss_process(distance=0)     # locate the duplication points on the coalescent tree
    print('\ngene_tree events:')
    pprint.pprint(events)
    print('\ngene_tree duplication_subtree:')
    qgtree.duplication_subtree(coalescent_process=coalescent_process, events=events)        # generate duplication subtrees


    # print('\nHAHAHAHAHAHA')
    # print(qstree.node_by_id(10))
    # print(qstree.node_by_name('E'))
    # print(qgtree.node_by_id(10))
    # print(qgtree.node_by_name('3*'))
    # print(qstree.node_by_id(0))
    # print(qgtree.node_by_id(14))

    # qgtree = qiuyi_tree.GeneTree(time_sequences=time_sequences, species_tree=qstree)        # construct newick coalescent tree
    # print('\ngene_tree ascii_art:')
    # print(qgtree.skbio_tree.ascii_art())
    # print('\ngene_nodes:')
    # qgtree.print_nodes()
    # print('\ngene_tree dup_loss_process:')
    # dup_events = qgtree.dup_loss_process(lambda_dup=0.2, lambda_loss=0.05, lambda_trans = 0.1)[0]      # locate the duplication points on the coalescent tree
    # trans_events = qgtree.dup_loss_process(lambda_dup=0.2, lambda_loss=0.05, lambda_trans = 0.1)[1]
    # print('\ngene_tree dup_events:')
    # pprint.pprint(dup_events)
    # print('\ngene_tree trans_events:')
    # pprint.pprint(trans_events)
    # print('\ngene_tree duplication_subtree:')
    # qgtree.duplication_subtree(coalescent_process=coalescent_process, dup_events=dup_events)        # generate duplication subtrees

    return
    

if __name__ == "__main__":
    main()
