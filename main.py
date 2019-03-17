import qiuyi_tree
import pprint


def main():
    qstree = qiuyi_tree.SpeciesTree(newick_path='data/tree_sample1.txt')
    print('\nsecies_tree ascii_art:')
    print(qstree.skbio_tree.ascii_art())
    print('\nspecies_nodes:')
    qstree.print_nodes()
    print('\ncoalescent:')
    coalescent_process = qstree.coalescent(distance_at_root=10000, lambda0=0.3)
    print('\ncoalescent_process:')
    pprint.pprint(coalescent_process)
    print('\ntime_sequences:')
    time_sequences = qstree.time_sequences(coalescent_process=coalescent_process)
    pprint.pprint(time_sequences)

    qgtree = qiuyi_tree.GeneTree(time_sequences=time_sequences, species_tree=qstree)
    print('\ngene_tree ascii_art:')
    print(qgtree.skbio_tree.ascii_art())
    print('\ngene_nodes:')
    qgtree.print_nodes()
    print('\ngene_tree dup_loss_process:')
    dup_events = qgtree.dup_loss_process(lambda_dup=0.1, lambda_loss=0.03)
    print('\ngene_tree dup_events:')
    pprint.pprint(dup_events)
    print('\ngene_tree duplication_subtree:')
    qgtree.duplication_subtree(coalescent_process=coalescent_process, dup_events=dup_events)


    # print('\nsubtree coalescent:')
    # sub_nodes = qstree.nodes.copy()
    # del sub_nodes[3]
    # del sub_nodes[5]
    # qstree_2 = qiuyi_tree.SpeciesTree(nodes=sub_nodes)
    # coalescent_process_2 = qstree_2.coalescent(distance_at_root=1.7, lambda0=1)
    # print('\nsubtree coalescent_process:')
    # pprint.pprint(coalescent_process_2)
    # print('\nsubtree time_sequences:')
    # time_sequences_2 = qstree_2.time_sequences(coalescent_process=coalescent_process_2)
    # pprint.pprint(time_sequences_2)

    return
    

if __name__ == "__main__":
    main()
