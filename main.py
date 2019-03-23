import qiuyi_tree
import pprint


def main():
    qstree = qiuyi_tree.SpeciesTree(newick_path='data/tree_sample.txt')    # read newick species tree
    print('\nsecies_tree ascii_art:')
    print(qstree.skbio_tree.ascii_art())
    print('\nspecies_nodes:')
    qstree.print_nodes()
    print('\ncoalescent:')
    coalescent_process = qstree.coalescent(distance_above_root=10000, lambda0=0.3)      # do coalescece based on the species tree
    print('\ncoalescent_process:')
    pprint.pprint(coalescent_process)
    print('\ntime_sequences:')
    time_sequences = qstree.time_sequences(coalescent_process=coalescent_process)       # convert to time sequence structure
    pprint.pprint(time_sequences)
    
    # print('\nTEST')
    # pprint.pprint(qstree.filter_coal_process(coalescent_process, '1*2*'))

    qgtree = qiuyi_tree.GeneTree(time_sequences=time_sequences, species_tree=qstree)        # construct newick coalescent tree
    print('\ngene_tree ascii_art:')
    print(qgtree.skbio_tree.ascii_art())
    print('\ngene_nodes:')
    qgtree.print_nodes()
    print('\ngene_tree DLT_process:')
    events = qgtree.dup_loss_process(lambda_dup=0.2, lambda_loss=0.2, lambda_trans=0.2)     # locate the duplication points on the coalescent tree
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
