import qiuyi_tree

def main():
    qstree = qiuyi_tree.SpeciesTree(lambda0=0.3, newick_path='data/tree_sample.txt')

    print('\nsecies_tree ascii_art:')
    print(qstree.skbio_tree.ascii_art())
    print('\nspecies_nodes:')
    qstree.print_nodes()
    print('\ncoalescent:')
    qstree.coalescent()

    qgtree = qiuyi_tree.GeneTree(species_tree=qstree)

    print('\ngene_tree ascii_art:')
    print(qgtree.skbio_tree.ascii_art())
    print('\ngene_nodes:')
    qgtree.print_nodes()

    return
    

if __name__ == "__main__":
    main()
