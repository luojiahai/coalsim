import qiuyi_tree

def main():
    qstree = qiuyi_tree.SpeciesTree(newick_path='data/tree_sample.txt',
                                    lambda0=0.3)

    print('\nsecies_tree ascii_art:')
    print(qstree.skbio_tree.ascii_art())
    print('\nspecies_nodes:')
    qstree.print_nodes()
    print('\ncoalescent:')
    qstree.coalescent()

    qgtree = qiuyi_tree.GeneTree(species_tree=qstree,
                                 lambda_dup=1,
                                 lambda_loss=0.3)

    print('\ngene_tree ascii_art:')
    print(qgtree.skbio_tree.ascii_art())
    print('\ngene_nodes:')
    qgtree.print_nodes()
    print('\ngene_tree dl_process:')
    qgtree.dl_process()

    return
    

if __name__ == "__main__":
    main()
