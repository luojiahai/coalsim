from shared_utility import *
from species_tree import *
from gene_tree import *

import pprint
import os, shutil
import numpy as np
import skbio
import time
import getopt, sys


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
        subtree = skbio.read(f, format='newick', into=skbio.tree.TreeNode)
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
    Debug.log_file = open('./output/log.txt', 'w')
    Debug.log(header='Log created on ' + time.ctime() + '\n')

    qstree = SpeciesTree(newick_path='data/tree_sample1.txt')    # read newick species tree
    Debug.save_tree_nodes(nodes=qstree.nodes, path='output/species_nodes_table.txt')

    SpeciesTree.global_species_tree = qstree
    SpeciesTree.lambda_coal = np.random.gamma(shape=3, scale=0.1, size=len(qstree.leaves))

    Debug.log(header='\nspecies_tree ascii_art:\n',
                         bodies=[qstree.skbio_tree.ascii_art()])
    Debug.log(header='\nspecies_nodes:\n',
                         bodies=qstree.nodes)
    Debug.log(header='\ncoalescent:\n')
    coalescent_process, _ = qstree.coalescent(distance_above_root=10000)      # do coalescece based on the species tree
    Debug.log(header='\ncoalescent_process:\n',
                         bodies=[coalescent_process],
                         pformat=True)
    time_sequences = qstree.time_sequences(coalescent_process=coalescent_process)       # convert to time sequence structure
    Debug.log(header='\ntime_sequences:\n',
                         bodies=[time_sequences],
                         pformat=True)
    

    qgtree = GeneTree(time_sequences=time_sequences, species_tree=qstree, coalescent_process=coalescent_process)        # construct newick coalescent tree

    GeneTree.lambda_dup = np.random.gamma(shape=1, scale=0.1, size=len(qgtree.leaves))
    GeneTree.lambda_loss = np.random.gamma(shape=1, scale=0.1, size=len(qgtree.leaves))
    GeneTree.lambda_trans = np.random.gamma(shape=1, scale=0.1, size=len(qgtree.leaves))
    GeneTree.recombination = options['recombination']
    GeneTree.hemiplasy = options['hemiplasy']

    Debug.save_tree_nodes(nodes=qgtree.nodes, path='output/gene_nodes_table.txt')
    # Debug.save_output(contents=[qgtree.skbio_tree],
    #                              path='output/newick_gene_subtrees/gene_tree.txt')
    Debug.log(header='\ngene_tree ascii_art:\n',
                         bodies=[qgtree.skbio_tree.ascii_art()])
    Debug.log(header='\ngene_nodes:\n',
                         bodies=qgtree.nodes)
    Debug.log(header='\ngene_tree dlt_process:\n')
    events = qgtree.dlt_process(distance=0)     # locate the duplication points on the coalescent tree
    Debug.log(header='\ngene_tree events:\n',
                         bodies=[events],
                         pformat=True)
    Debug.log(header='\ngene_tree dt_subtree:\n')
    qgtree.dt_subtree(coalescent_process=coalescent_process, events=events, path=tree_path)        # generate duplication subtrees


    final_result = qgtree.skbio_tree.deepcopy()
    loss_nodes = build_tree(final_result, './output/tree')
    Debug.save_output(contents=[final_result,final_result.ascii_art()],
                                 path='./output/final_result.txt')
    final_result_cut = cut_tree(final_result, loss_nodes)
    Debug.save_output(contents=[final_result_cut,final_result_cut.ascii_art()],
                                 path='./output/final_result_cut.txt')
    Debug.log_file.close()
    
    Debug.save_output(contents=[final_result_cut],path='./output/final_result_cut_newick.txt')

    Debug.save_output(contents='species_tree: \n',path='./output/final_result_details.txt')
    qgtree.construct_final_gene_nodes(final_result_cut)
    # Debug.save_output(contents=qgtree.nodes,path='./output/final_result_details.txt')
        # print(node)
    # gene_root = qgtree.nodes[-1]
    qgtree.post_order_fake_id()
    l = []
    for node in qgtree.nodes:
        l.append(str('real_id: '+ str(node.node_id) + ' fake_id: '+ str(node.fake_node_id)))
    Debug.save_output(contents=qgtree.nodes + l, path='./output/final_result_details.txt')


    # species_root = SpeciesTree.global_species_tree.root
    # SpeciesTree.global_species_tree.post_order_fake_id(species_root)
    # for node in SpeciesTree.global_species_tree.nodes:
    #     print('real_id: ', node.node_id, ' fake_id: ', node.fake_node_id)


    print('Number of events: ')
    print('Duplicaton: ' + str(Debug.event_count['d']))
    print('Loss: ' + str(Debug.event_count['l']))
    print('Transfer: ' + str(Debug.event_count['t']))
    print('ILS: ' + str(Debug.event_count['i']))

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


def parse_arg(argv):
    try: 
        opts, args = getopt.getopt(argv,'r:h:',
        ['recombination=','hemiplasy=','help'])
    except getopt.GetoptError:
        print('Usage: {} -r <recombination> -h <hemiplasy>'.format(sys.argv[0]))
        sys.exit()
    if(opts):
        for opt, arg in opts:
            if opt in ('-r', '--recombination'):
                recombination = arg
            elif opt in ('-h', '--hemiplasy'):
                hemiplasy = arg
            elif opt in ('--help'):
                print('Usage: {} -r <recombination> -h <hemiplasy>'.format(sys.argv[0]))
                sys.exit()
    else:
        recombination = 1
        hemiplasy = 1

    return recombination, hemiplasy


if __name__ == '__main__':
    recombination, hemiplasy = parse_arg(sys.argv[1:])
    options = {
        'recombination': int(recombination),
        'hemiplasy': int(hemiplasy)
    }
    main(options)
