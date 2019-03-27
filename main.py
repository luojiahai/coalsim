import qiuyi_tree
import pprint
import os, shutil
import numpy as np
import skbio


def dir_recurse(gene_tree, path):
    files = os.listdir(path)
    for f in files:
        f_path = os.path.join(path, f)
        if os.path.isdir(f_path):
            dir_recurse(gene_tree, f_path)
            
    _id = '_' + path.split('_')[-1]

    subtree_path = os.path.join(path, 'gene_tree.txt')
    f = open(subtree_path)
    subtree = skbio.read(f, format="newick", into=skbio.tree.TreeNode)
    f.close()

    event_path = os.path.join(path, 'event.txt')
    f = open(event_path)
    line = f.readline()
    splited = line.strip().split(',')
    name = splited[0]
    distance = float(splited[1])

    new_node = skbio.TreeNode()
    child = gene_tree.find(name)
    parent = child.parent
    new_node.name = name + _id
    new_node.length = child.length - distance
    new_node.parent = parent
    new_node.children.append(child)
    child.length = distance
    child.parent = new_node
    for i in range(len(parent.children)):
        if (parent.children[i].name == child.name):
            del parent.children[i]
            break
    parent.children.append(new_node)
    new_node.children.append(subtree)
    subtree.parent = new_node
    for node in subtree.traverse():
        node.name = node.name + _id

def build_tree(gene_tree, path):
    files = os.listdir(path)
    for f in files:
        f_path = os.path.join(path, f)
        if os.path.isdir(f_path):
            dir_recurse(gene_tree, f_path)
    return

def main():
    shutil.rmtree('./output')
    os.mkdir('./output')
    os.mkdir('./output/newick_gene_subtrees')
    os.mkdir('./output/subtrees')
    tree_path = './output/tree'
    os.mkdir(tree_path)

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
    qiuyi_tree.GeneTree.lambda_loss = np.random.gamma(shape=1, scale=0.1, size=len(qgtree.leaves))
    qiuyi_tree.GeneTree.lambda_trans = np.random.gamma(shape=1, scale=0.1, size=len(qgtree.leaves))

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
    qgtree.duplication_subtree(coalescent_process=coalescent_process, events=events, path=tree_path)        # generate duplication subtrees


    final = qgtree.skbio_tree.deepcopy()
    build_tree(final, './output/tree')
    f = open('./output/final.txt', 'w')
    f.write(str(final))
    f.write(str(final.ascii_art()))
    f.close

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
