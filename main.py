import qiuyi_tree

def main():
    qiuyi_tree.newick_to_table(input_path='data/tree_sample1.txt', output_path='data/nodes_table.txt')

    qstree = qiuyi_tree.SpeciesTree(table_file_path='data/nodes_table.txt', lambda0=0.3)
    qstree.print_nodes()
    qstree.coalescent()

    qgtree = qiuyi_tree.GeneTree(species_tree=qstree)
    

if __name__ == "__main__":
    main()
