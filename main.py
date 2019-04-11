import qiuyi_tree
import pprint
import os, shutil
import numpy as np
import skbio
import time


def build_tree_recurse(gene_tree, path):
    loss_nodes = []
    current_tree = gene_tree
    path_splited = path.split('_')
    _id = '_' + path.split('_')[-1] if len(path_splited) > 1 else ''
    parent_name_splited = gene_tree.name.split('_')
    _parent_id = '_' + parent_name_splited[-1] if len(parent_name_splited) > 1 else ''

    if (_id):
        subtree_path = os.path.join(path, 'gene_tree.txt')
        f = open(subtree_path)
        subtree = skbio.read(f, format="newick", into=skbio.tree.TreeNode)
        f.close()
        current_tree = subtree

        event_path = os.path.join(path, 'event.txt')
        f = open(event_path)
        line = f.readline()
        splited = line.strip().split(',')
        node_name = splited[0]
        distance = float(splited[1])
        event_name = '_' + splited[2][0]
        f.close()

        new_dt_node = skbio.TreeNode()
        child = None
        for node in gene_tree.traverse():
            if node.name == (node_name + _parent_id):
                child = node
        parent = child.parent
        new_dt_node.name = node_name + event_name + _id
        new_dt_node.length = child.length - distance
        new_dt_node.parent = parent
        new_dt_node.children.append(child)
        child.length = distance
        child.parent = new_dt_node
        for i in range(len(parent.children)):
            if (parent.children[i].name == child.name):
                del parent.children[i]
                break
        parent.children.append(new_dt_node)
        new_dt_node.children.append(subtree)
        subtree.parent = new_dt_node
        for node in subtree.traverse():
            node.name = node.name + _id

    files = os.listdir(path)
    for f in files:
        file_path = os.path.join(path, f)
        if os.path.isdir(file_path):
            loss_nodes += build_tree_recurse(current_tree, file_path)

    files = os.listdir(path)
    for f in files:
        file_name = f.split('_')
        if (file_name and file_name[0] == 'loss'):
            loss_path = os.path.join(path, f)
            file_ = open(loss_path)
            line = file_.readline()
            splited = line.split(',')
            node_l_name = splited[0]
            node_l_distance = float(splited[1])

            new_l_node = skbio.TreeNode()
            child = None
            for node in current_tree.traverse():
                if node.name == (node_l_name + _id):
                    child = node
            parent = child.parent
            new_l_node.name = node_l_name + '_l' + _id
            new_l_node.length = child.length - node_l_distance
            new_l_node.parent = parent
            new_l_node.children.append(child)
            child.length = node_l_distance
            child.parent = new_l_node
            for i in range(len(parent.children)):
                if (parent.children[i].name == child.name):
                    del parent.children[i]
                    break
            parent.children.append(new_l_node)
            loss_nodes.append(new_l_node)

            file_.close()

    for f in files:
        file_name = f.split('_')
        if (file_name and file_name[0] == 'ils'):
            # do something
            ils_path = os.path.join(path, f)
            file_ = open(ils_path)
            line = file_.readline()
            splited = line.split(',')
            node_name = splited[0]
            for node in current_tree.traverse():
                if node.name == (node_name + _id):
                    node.name += '_i'
                    break
            file_.close()

    return loss_nodes

def build_tree(gene_tree, path):
    return build_tree_recurse(gene_tree, path)

def cut_tree(final_tree, loss_nodes):
    final = final_tree.deepcopy()
    for node in loss_nodes:
        final.remove_deleted(lambda x: x.name == node.name)
    final.prune()
    return final

def main(options):
    shutil.rmtree('./output')
    os.mkdir('./output')
    # os.mkdir('./output/newick_gene_subtrees')
    # os.mkdir('./output/subtrees')
    tree_path = './output/tree'
    os.mkdir(tree_path)
    qiuyi_tree.Debug.log_file = open('./output/log.txt', 'w')
    qiuyi_tree.Debug.log(header='Log created on ' + time.ctime() + '\n')

    qstree = qiuyi_tree.SpeciesTree(newick_path='data/tree_sample1.txt')    # read newick species tree
    qiuyi_tree.Debug.save_tree_nodes(nodes=qstree.nodes, path='output/species_nodes_table.txt')

    qiuyi_tree.SpeciesTree.global_species_tree = qstree
    qiuyi_tree.SpeciesTree.lambda_coal = qiuyi_tree.np.random.gamma(shape=3, scale=0.1, size=len(qstree.leaves))

    qiuyi_tree.Debug.log(header='\nspecies_tree ascii_art:\n',
                         bodies=[qstree.skbio_tree.ascii_art()])
    qiuyi_tree.Debug.log(header='\nspecies_nodes:\n',
                         bodies=qstree.nodes)
    qiuyi_tree.Debug.log(header='\ncoalescent:\n')
    coalescent_process, _ = qstree.coalescent(distance_above_root=10000)      # do coalescece based on the species tree
    qiuyi_tree.Debug.log(header='\ncoalescent_process:\n',
                         bodies=[coalescent_process],
                         pformat=True)
    time_sequences = qstree.time_sequences(coalescent_process=coalescent_process)       # convert to time sequence structure
    qiuyi_tree.Debug.log(header='\ntime_sequences:\n',
                         bodies=[time_sequences],
                         pformat=True)
    

    qgtree = qiuyi_tree.GeneTree(time_sequences=time_sequences, species_tree=qstree)        # construct newick coalescent tree

    qiuyi_tree.GeneTree.lambda_dup = qiuyi_tree.np.random.gamma(shape=1, scale=0.1, size=len(qgtree.leaves))
    qiuyi_tree.GeneTree.lambda_loss = qiuyi_tree.np.random.gamma(shape=1, scale=0.1, size=len(qgtree.leaves))
    qiuyi_tree.GeneTree.lambda_trans = qiuyi_tree.np.random.gamma(shape=1, scale=0.1, size=len(qgtree.leaves))
    qiuyi_tree.GeneTree.dup_recombination = options['dup_recombination']
    qiuyi_tree.GeneTree.trans_hemiplasy = options['trans_hemiplasy']

    qiuyi_tree.Debug.save_tree_nodes(nodes=qgtree.nodes, path='output/gene_nodes_table.txt')
    # qiuyi_tree.Debug.save_output(contents=[qgtree.skbio_tree],
    #                              path='output/newick_gene_subtrees/gene_tree.txt')
    qiuyi_tree.Debug.log(header='\ngene_tree ascii_art:\n',
                         bodies=[qgtree.skbio_tree.ascii_art()])
    qiuyi_tree.Debug.log(header='\ngene_nodes:\n',
                         bodies=qgtree.nodes)
    qiuyi_tree.Debug.log(header='\ngene_tree dlt_process:\n')
    events = qgtree.dlt_process(distance=0)     # locate the duplication points on the coalescent tree
    qiuyi_tree.Debug.log(header='\ngene_tree events:\n',
                         bodies=[events],
                         pformat=True)
    qiuyi_tree.Debug.log(header='\ngene_tree dt_subtree:\n')
    qgtree.dt_subtree(coalescent_process=coalescent_process, events=events, path=tree_path)        # generate duplication subtrees


    final_result = qgtree.skbio_tree.deepcopy()
    loss_nodes = build_tree(final_result, './output/tree')
    qiuyi_tree.Debug.save_output(contents=[final_result,final_result.ascii_art()],
                                 path='./output/final_result.txt')
    final_result_cut = cut_tree(final_result, loss_nodes)
    qiuyi_tree.Debug.save_output(contents=[final_result_cut,final_result_cut.ascii_art()],
                                 path='./output/final_result_cut.txt')

    qiuyi_tree.Debug.log_file.close()

    print('Number of events: ')
    print('Duplicaton: ' + str(qiuyi_tree.Debug.event_count['d']))
    print('Loss: ' + str(qiuyi_tree.Debug.event_count['l']))
    print('Transfer: ' + str(qiuyi_tree.Debug.event_count['t']))
    print('ILS: ' + str(qiuyi_tree.Debug.event_count['i']))

    # for node in final_result_cut.tips():
    #     print(final_result_cut.distance(node))

    # print('\nTest')
    # print(qstree.node_by_id(10))
    # print(qstree.node_by_name('E'))
    # print(qgtree.node_by_id(10))
    # print(qgtree.node_by_name('3*'))
    # print(qstree.node_by_id(0))
    # print(qgtree.node_by_id(14))
    return
    

if __name__ == "__main__":
    options = {
        'dup_recombination': 0,
        'trans_hemiplasy': 0
    }
    main(options)
