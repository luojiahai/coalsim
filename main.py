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
        event_index = '_' + splited[3]
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
                splited = node.name.split('_')
                if (_id == ''):
                    if (splited[0] == node_l_name):
                        child = node
                else:
                    if (splited[0] == node_l_name and ('_' + splited[-1]) == _id):
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
            _index = '_' + file_name[1]
            ils_path = os.path.join(path, f)
            file_ = open(ils_path)
            line = file_.readline()
            splited = line.split(',')
            node_name = splited[0]
            for node in current_tree.traverse():
                if node.name == (node_name + _id):
                    node.name += '_i' + _index
                    break
            file_.close()

    for f in files:
        file_name = f.split('_')
        if (file_name and file_name[0] == 's'):
            _index = '_' + file_name[1]
            s_path = os.path.join(path, f)
            file_ = open(s_path)
            line = file_.readline()
            splited = line.split(',')
            node_name = splited[0]
            for node in current_tree.traverse():
                if node.name == (node_name + _id):
                    node.name += '_s' + _index
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
    Debug.summary_file = open('./output/summary.txt', 'w')
    Debug.summary(header='Summary created on ' + time.ctime() + '\n')
    Debug.summary(header='Events:\n') 
    Debug.summary(header='type\tspecies_node\tclade_set\n')   

    qstree = SpeciesTree(newick_path='data/tree_sample1.txt')    # read newick species tree
    Debug.save_tree_nodes(nodes=qstree.nodes, path='output/species_nodes_table.txt')

    SpeciesTree.global_species_tree = qstree
    # I want fake id
    SpeciesTree.global_species_tree.post_order_fake_id()
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

    GeneTree.lambda_dup = np.random.gamma(shape=0.5, scale=0.1, size=len(qgtree.leaves))
    GeneTree.lambda_loss = np.random.gamma(shape=0.5, scale=0.1, size=len(qgtree.leaves))
    GeneTree.lambda_trans = np.random.gamma(shape=0.3, scale=0.1, size=len(qgtree.leaves))
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
    
    Debug.log(header='\nfull_events:\n',
                         bodies=qgtree.full_events,
                         pformat=True)
    # save final gene tree
    final_result = qgtree.skbio_tree.deepcopy()
    loss_nodes = build_tree(final_result, './output/tree')
    # for node in loss_nodes:
    #     if (node.parent.parent):
    #         node.parent.parent.name += '_lp'
    #     else: 
    #         node.parent.event = 'l'
    Debug.save_output(contents=[final_result,final_result.ascii_art()],
                                 path='./output/final_result.txt')
    final_result_cut = cut_tree(final_result, loss_nodes)
    Debug.save_output(contents=[final_result_cut,final_result_cut.ascii_art()],
                                 path='./output/final_result_cut.txt')
    Debug.save_output(contents=[final_result_cut],path='./output/final_result_cut_newick.txt')

    if (not final_result_cut):
        print("EXCEPTION: ALL LOST")
        Debug.log_file.close()
        Debug.summary_file.close()
        return
            
    # species tree table
    # SpeciesTree.global_species_tree.post_order_fake_id()

    Debug.summary(header='\nSpecies_tree_table:\n')
    Debug.summary(header='node_id\tclade_set\n')

    for node in SpeciesTree.global_species_tree.nodes:
        clade = node.clade
        for i in range(len(clade)):
            clade[i] = SpeciesTree.global_species_tree.get_fake_id_from_real_id(clade[i])
        Debug.summary(header=str(node.node_id) + '\t'
                                + str(node.fake_node_id) + '\t'
                                + str(clade) + '\n')

    if (not final_result_cut):
        print("EXCEPTION: ALL LOST")
    else:
        # gene tree table
        qgtree.construct_final_gene_nodes(final_result_cut)
        qgtree.post_order_fake_id()
        
        Debug.summary(header='\nGene_tree_table:\n')
        Debug.summary(header='gene_node_id\tclade_set\tevent\tspecies_node_id\n')
        for node in qgtree.nodes:      
            clade = node.name.split('*')[:-1]
            for i in range(len(clade)):
                clade[i] = SpeciesTree.global_species_tree.get_fake_id_from_real_id(clade[i])                     
            Debug.summary(header=str(node.node_id) + '\t'
                                + str(node.fake_node_id) + '\t'
                                + str(clade) + '\n' )
    
        Debug.summary(header='\nevent_table:\n')
        for node in qgtree.nodes:
            if ('_d' in node.name):
                node.event = 'd'
            elif ('_t' in node.name):
                node.event = 't'
            # elif ('_l' in node.name):
            #     node.event = 'l'
            elif ('_i' in node.name):
                node.event = 'i'
            elif ('_s' in node.name):
                node.event = 's'
            index = node.name.split('_')[-1]
            find_it = False
            if (node.event):
                for event in qgtree.full_events:
                    if (str(event['index'])==index):
                        species_fake_id = SpeciesTree.global_species_tree.get_fake_id_from_real_id(event['species_node_id'])
                        find_it = True
                        break
            if (find_it == False):
                species_fake_id = -1
            clade = node.name.split('*')[:-1]
            for i in range(len(clade)):
                clade[i] = SpeciesTree.global_species_tree.get_fake_id_from_real_id(clade[i])
            if (node.event):
                Debug.event_count[node.event] += 1
                Debug.summary(header=str(node.node_id) + '\t'
                                    + str(node.fake_node_id) + '\t'
                                    + str(clade) + '\t' 
                                    + str(node.event) + '\t'
                                    + str(species_fake_id) + '\n' )

    Debug.summary(header='\ntime_sequences:\n',
                         bodies=qgtree.full_events,
                         pformat=True)

    print('Number of events: ')
    print('Speciation: ' + str(Debug.event_count['s']))
    print('Duplicaton: ' + str(Debug.event_count['d']))
    print('Loss: ' + str(Debug.event_count['l']))
    print('Transfer: ' + str(Debug.event_count['t']))
    print('ILS: ' + str(Debug.event_count['i']))

    # node = final_result_cut.tips()[0]
    # distance_to_root = final_result_cut.distance(node)
    # for node in final_result_cut.tips():
    #     if (distance_to_root != final_result_cut.distance(node)):
    #         print("invalid output!")
    #         break
    #         return

    Debug.log_file.close()
    Debug.summary_file.close()
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
