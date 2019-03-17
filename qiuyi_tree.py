from io import StringIO
import skbio
import skbio.tree

import numpy as np
import random
import collections

import pprint

np.random.seed(4)


class TreeNode(object):
    def __init__(self,
                 node_id=None,
                 name=None,
                 parent=None,
                 distance_to_parent=None,
                 children=None):
        self.node_id = node_id
        self.name = name
        self.parent = parent
        self.distance_to_parent = distance_to_parent
        self.children = children if children else []
        self.distance_to_children = []
        return
    
    def __repr__(self):
        return "node_id: {}, name: {}, parent: {}, distance_to_parent: {}, children: {}, distance_to_children: {}".format(
                self.node_id, 
                self.name,
                self.parent, 
                self.distance_to_parent,
                self.children,
                self.distance_to_children)


class GenericTree(object):
    def __init__(self):
        self.skbio_tree = None
        self.nodes = []
        self.nodes_name_dict = {}
        self.root = None
        self.leaves = []
        return

    def save_to_file(self, path):
        f = open(path, 'w')
        f.write('id\tname\tparent\td2p\n')
        for node in self.nodes:
            f.write(str(node.node_id) + '\t' + node.name + '\t' + str(node.parent) + '\t' + str(node.distance_to_parent) + '\n')
        f.close()

    def newick_to_table(self,
                        output_path, 
                        input_path=None, 
                        skbio_tree=None):
        tree = None
        if (skbio_tree):
            tree = skbio_tree
        else:
            f = open(input_path)
            tree = skbio.read(f, format="newick", into=skbio.tree.TreeNode)
            f.close()

        def _parse(tree):
            node = {
                'object': tree,
                'name': tree.name,
                'parent': tree.parent,
                'children': [],
                'distance': tree.length
            }
            if tree.is_tip():
                return node
            for children in tree.children:
                node['children'].append(_parse(children))
            return node

        def _rename(node):
            ret = node.copy()
            name = ''
            for i in range(len(ret['children'])):
                ret['children'][i] = _rename(ret['children'][i])
                name += ret['children'][i]['name']
            if (not ret['name']):
                ret['name'] = name
            return ret

        def _to_list(node, root):
            d = node.copy()
            ret = []
            for i in range(len(d['children'])):
                ret += _to_list(d['children'][i], root=root)
            del d['children']
            d['distance_to_root'] = d['object'].distance(root)
            ret.append(d)
            if (d['object'] is root):
                ret = sorted(ret, key=lambda x: x['distance_to_root'], reverse=True)
            return ret

        root = _parse(tree)
        root = _rename(root)
        nodes = _to_list(root, root=root['object'])

        # output to file
        f = open(output_path, 'w')
        f.write('id\tname\tparent\td2p\n')
        for i in range(len(nodes)):
            parent_index = 'None'
            for j in range(len(nodes)):
                if (nodes[i]['parent'] is nodes[j]['object']):
                    parent_index = str(j)
            f.write(str(i) + '\t' + nodes[i]['name'] + '\t' + parent_index + '\t' + str(nodes[i]['distance']) + '\n')
        f.close()

        return tree

    def print_nodes(self):
        for node in self.nodes:
            print(node)
        return

    def __process_tree_recurse(self, tree):
        if (tree.name):
            return tree.name
        else:
            name = ''
            for child in tree.children:
                name += self.__process_tree_recurse(child)
            tree.name = name
            return tree.name

    def construct_nodes(self, path, process_tree=False):
        # open table file and construct tree nodes
        f = open(path, 'r')
        f.readline()
        for line in f:
            splited = line.strip().split('\t')
            node_id = int(splited[0])
            name = splited[1]
            parent = splited[2]
            if (parent == 'None'):
                parent = -1
            d2p = splited[3]
            if (d2p == 'None'):
                d2p = -1.0
            self.nodes.append(TreeNode(node_id=node_id,
                                       name=name,
                                       parent=int(parent),
                                       distance_to_parent=float(d2p)))
        # process tree
        if (process_tree):
            self.__process_tree_recurse(self.skbio_tree)
        # create dict
        for node in self.nodes:
            self.nodes_name_dict[node.name] = node
        # find children
        for node in self.nodes:
            children, node.distance_to_children = self.children_distances(self.skbio_tree, node.name)
            for child in children:
                node_id = self.nodes_name_dict[child.name].node_id
                node.children.append(node_id)
        return

    def __children_distances_recurse(self, tree, name):
        ret = None
        if (tree.name == name):
            distances = []
            for child in tree.children:
                distances.append(tree.distance(child))
            ret = (tree.children, distances)
            return ret
        else:
            for child in tree.children:
                ret = self.__children_distances_recurse(child, name)
                if (ret): return ret

    def children_distances(self, tree, name):
        return self.__children_distances_recurse(tree, name)
        

class SpeciesTree(GenericTree):
    def __init__(self,
                 newick_path=None,
                 nodes=None):
        GenericTree.__init__(self)
        if (not newick_path):
            self.nodes = nodes
        else:
            self.__construct_species_nodes(newick_path)

        max_node_id = -1
        for node in self.nodes:
            if (node.node_id > max_node_id):
                max_node_id = node.node_id
                self.root = node
                
        self.nodes_id_dict = {}
        for node in self.nodes:
            self.nodes_id_dict[node.node_id] = node
        
        self.leaves = [node.node_id for node in self.nodes if not node.children]
        self.total_distance = self.__distance_to_root_recurse(node_id=self.leaves[0])
        return

    def __construct_species_nodes(self, newick_path):
        output_path = 'output/species_nodes_table.txt'
        self.skbio_tree = super().newick_to_table(input_path=newick_path, output_path=output_path)
        super().construct_nodes(output_path, process_tree=True)
        return

    def __star_sorted(self, couple):
        string = ''
        for e in couple:
            string += e
        splited = string.split('*')[:-1]
        splited = sorted([int(e) for e in splited])
        return [str(e) + '*' for e in splited]

    def __coalescent_recurse(self,
                             node_id,
                             distance,
                             clade_set,
                             lambda0,
                             coalescent_process):
        if (len(clade_set[node_id]) <= 1):
            return
        else:
            lambda_c = len(clade_set[node_id]) * lambda0
            distance_fake = np.random.exponential(scale=1.0/lambda_c)
            if (distance < distance_fake):
                return
            else:
                if (len(clade_set[node_id]) >= 2):
                    temp_set = sorted(clade_set[node_id])
                    couple = np.random.choice(clade_set[node_id], size=2, replace=False)
                    clade_set[node_id] = [''.join(self.__star_sorted(couple))] + [e for e in clade_set[node_id] if e not in couple]

                    # print process
                    print("initial node " + str(node_id) + ": " + str(temp_set))
                    print("coalescent at node " + str(node_id) + ": " + str(clade_set[node_id]) + ", " + "distance = " + str(distance_fake))

                    # save process
                    coalescent_process[str(node_id)].append({
                        'from_set': temp_set, 
                        'to_set': clade_set[node_id].copy(),
                        'distance': distance_fake
                    })
                else:
                    return
                distance = distance - distance_fake
                self.__coalescent_recurse(node_id=node_id, 
                                          distance=distance, 
                                          clade_set=clade_set, 
                                          lambda0=lambda0,
                                          coalescent_process=coalescent_process)
        return

    def coalescent(self, distance_at_root, lambda0):
        nodes = self.nodes
        root = self.root
        coalescent_process = collections.defaultdict(list)

        old_leaves = [node.node_id for node in nodes if not node.children]
        new_leaves = []
        clade_set = {}
        labelled = {}
        nodes_id_dict = {}
        for node in nodes:
            labelled[node.node_id] = False
            nodes_id_dict[node.node_id] = node
            clade_set[node.node_id] = [str(node.node_id) + '*'] if not node.children else []

        while (True):
            for leaf in old_leaves:
                if (leaf == root.node_id):
                    self.__coalescent_recurse(node_id=root.node_id, 
                                              distance=distance_at_root, 
                                              clade_set=clade_set,
                                              lambda0=lambda0,
                                              coalescent_process=coalescent_process)
                    break
                else:
                    parent = nodes_id_dict[leaf].parent
                    children = nodes_id_dict[parent].children
                    if (labelled[leaf]):
                        continue
                    labelled[leaf] = True
                    if (len(clade_set[children[0]]) != 0 
                        and len(clade_set[children[1]]) != 0):
                        self.__coalescent_recurse(node_id=children[0], 
                                                  distance=nodes_id_dict[children[0]].distance_to_parent,
                                                  clade_set=clade_set,
                                                  lambda0=lambda0,
                                                  coalescent_process=coalescent_process)
                        self.__coalescent_recurse(node_id=children[1], 
                                                  distance=nodes_id_dict[children[1]].distance_to_parent,
                                                  clade_set=clade_set,
                                                  lambda0=lambda0,
                                                  coalescent_process=coalescent_process)
                        clade_set[parent] = list(set().union(clade_set[children[0]], clade_set[children[1]]))
                        if (len(new_leaves) > 0):
                            new_leaves = [e for e in new_leaves if e != children[0] and e != children[1]]
                        new_leaves.append(parent)
                    else:
                        new_leaves.append(leaf)
            if (leaf == root.node_id):
                break
            temp_new_leaves = []
            for new_leaf in new_leaves:
                if (new_leaf not in temp_new_leaves):
                    temp_new_leaves.append(new_leaf)
            old_leaves = temp_new_leaves.copy()
            new_leaves = []
            labelled = {}
            for node in nodes:
                labelled[node.node_id] = False
        return coalescent_process

    def __distance_to_root_recurse(self, node_id):
        if (self.nodes_id_dict[node_id].parent < 0 or 
            node_id == self.root.node_id):
            return 0
        else:
            d2p = self.nodes_id_dict[node_id].distance_to_parent
            parent = self.nodes_id_dict[node_id].parent
            return d2p + self.__distance_to_root_recurse(parent)
        
    def __walking_distance(self, node_id, branch_distance):
        return branch_distance + (self.total_distance - self.__distance_to_root_recurse(node_id))

    def __star_in_set(self, target, clade):
        if (len(target) < len(clade)):
            splited_target = target.split('*')[:-1]
            splited_clade = clade.split('*')[:-1]
            return set(splited_target).issubset(set(splited_clade))
        else:
            return False
    
    def __reverse_time_order(self, target_star, coalescent_process):
        sequence = []
        for k, v in coalescent_process.items():
            branch_distance = 0.0
            for elem in v:
                branch_distance += elem['distance']
                if (target_star in elem['from_set'] and target_star not in elem['to_set']):
                    for e in elem['to_set']:
                        if self.__star_in_set(target_star, e):
                            couple = e.replace(target_star, '')
                            walking_distance = self.__walking_distance(int(k), branch_distance=branch_distance)
                            # pair = (couple, walking_distance)
                            pair = (e, walking_distance)
                            sequence.append(pair)
                            sequence += self.__reverse_time_order(target_star=e, 
                                                                  coalescent_process=coalescent_process)
        return sequence

    def time_sequences(self, coalescent_process):
        time_sequences = {}
        for leaf in self.leaves:
            time_sequences[str(leaf)] = self.__reverse_time_order(target_star=str(leaf)+'*', 
                                                                  coalescent_process=coalescent_process)
        return time_sequences


class GeneTree(GenericTree):
    def __init__(self,
                 time_sequences,
                 species_tree):
        GenericTree.__init__(self)
        self.time_sequences = time_sequences
        self.__construct_gene_nodes()

        max_node_id = -1
        for node in self.nodes:
            if (node.node_id > max_node_id):
                max_node_id = node.node_id
                self.root = node
        
        self.species_tree = species_tree
        return

    def __star_replace(self, string, substring):
        a = string.split('*')[:-1]
        b = substring.split('*')[:-1]
        diff = set(a).difference(set(b))
        return ''.join([e + '*' for e in sorted(list(diff))])

    def __distance(self, node_name, parent_name):
        for leaf, sequence in self.time_sequences.items():
            if (len(node_name) == 2 and node_name[0] == leaf):
                for pair in sequence:
                    if (pair[0] == parent_name):
                        return pair[1]
            else:
                prev_pair = None
                for pair in sequence:
                    if (prev_pair != None and prev_pair[0] == node_name and pair[0] == parent_name):
                        return pair[1] - prev_pair[1]
                    prev_pair = pair
        return None

    def __create_skbio_tree_recurse(self, skbio_tree_node):
        # one node (leaf)
        if (len(skbio_tree_node.name) == 2):
            skbio_tree_node.length = self.__distance(skbio_tree_node.name, skbio_tree_node.parent.name)
            return
        # two nodes
        elif (len(skbio_tree_node.name) == 4):
            child_one_name = skbio_tree_node.name[:2]
            child_two_name = skbio_tree_node.name[2:]
            child_one = skbio.tree.TreeNode(name=child_one_name, 
                                            length=self.__distance(child_one_name, skbio_tree_node.name), 
                                            parent=skbio_tree_node)
            child_two = skbio.tree.TreeNode(name=child_two_name, 
                                            length=self.__distance(child_two_name, skbio_tree_node.name),
                                            parent=skbio_tree_node)
            skbio_tree_node.children = [child_one, child_two]
            return
        is_found = False
        for leaf, sequence in self.time_sequences.items():
            prev_pair = None
            for pair in sequence:
                if (prev_pair != None and skbio_tree_node.name == pair[0]):
                    child_one_name = prev_pair[0]
                    child_two_name = self.__star_replace(skbio_tree_node.name, prev_pair[0])
                    child_one = skbio.tree.TreeNode(name=child_one_name, 
                                                    length=self.__distance(child_one_name, skbio_tree_node.name), 
                                                    parent=skbio_tree_node)
                    child_two = skbio.tree.TreeNode(name=child_two_name, 
                                                    length=self.__distance(child_two_name, skbio_tree_node.name),
                                                    parent=skbio_tree_node)
                    self.__create_skbio_tree_recurse(child_one)
                    self.__create_skbio_tree_recurse(child_two)
                    skbio_tree_node.children = [child_one, child_two]
                    is_found = True
                    break
                prev_pair = pair
            if (is_found):
                break
        return

    def __construct_gene_nodes(self):
        # construct skbio tree from time sequence
        tree = skbio.tree.TreeNode()
        tree.name = next(iter(self.time_sequences.values()))[-1][0]
        self.__create_skbio_tree_recurse(tree)
        tree.length = None
        self.skbio_tree = tree
        super().newick_to_table(skbio_tree=tree, output_path='output/gene_nodes_table.txt')
        super().construct_nodes('output/gene_nodes_table.txt', process_tree=False)
        return

    def __dl_process_recurse(self, tree, distance, lambda_dup, lambda_loss, dup_events):
        node = self.nodes_name_dict[tree.name]
        distance_dup = np.random.exponential(scale=1.0/lambda_dup)
        distance_loss = np.random.exponential(scale=1.0/lambda_loss)
        if (distance_dup < distance_loss and distance_dup < distance):
            print('duplication at node ' + str(node.node_id) + ' (' + node.name + ')' + ' with distance ' + str(distance - distance_dup))
            dup_events.append({
                'node_id': node.node_id, 
                'name': node.name, 
                'distance': distance - distance_dup
            })
            self.__dl_process_recurse(tree, distance - distance_dup, lambda_dup, lambda_loss, dup_events)
        elif (distance_loss <= distance_dup and distance_loss < distance):
            print('loss at node ' + str(node.node_id) + ' (' + node.name + ')' + ' with distance ' + str(distance_loss))
        else:
            print('nothing happened at node ' + str(node.node_id) + ' (' + node.name + ')')
            if (node.children):
                child_one = tree.children[0]
                child_two = tree.children[1]
                distance_to_child_one = node.distance_to_children[0]
                distance_to_child_two = node.distance_to_children[1]
                self.__dl_process_recurse(child_one, distance_to_child_one, lambda_dup, lambda_loss, dup_events)
                self.__dl_process_recurse(child_two, distance_to_child_two, lambda_dup, lambda_loss, dup_events)
            else:
                print('reach the end of node ' + str(node.node_id) + ' (' + node.name + ')')
        return

    def dl_process(self, lambda_dup, lambda_loss):
        dup_events = []
        self.__dl_process_recurse(self.skbio_tree, 
                                  distance=0, 
                                  lambda_dup=lambda_dup, 
                                  lambda_loss=lambda_loss, 
                                  dup_events=dup_events)
        return dup_events
        
    def __duplication_subtree_recurse(self, event, node_id, coal_distance):
        species_skbio_tree = self.species_tree.skbio_tree
        name = self.species_tree.nodes_id_dict[node_id].name

        subtree = species_skbio_tree.find(name)
        subtree_names = [node.name for node in subtree.traverse()]
        subtree_nodes = [node for node in self.species_tree.nodes if node.name in subtree_names]

        species_subtree = SpeciesTree(nodes=subtree_nodes)
        species_subtree.skbio_tree = subtree
        print('\nspecies_subtree_nodes:')
        species_subtree.print_nodes()

        distance_at_root = event['distance'] + coal_distance
        print('\nspecies_subtree_coal:')
        # species_subtree_coal = species_subtree.coalescent(distance_at_root=distance_at_root, lambda0=0.3)
        species_subtree_coal_process = species_subtree.coalescent(distance_at_root=10000, lambda0=0.3)

        print('\nspecies_subtree_coal_process:')
        pprint.pprint(species_subtree_coal_process)

        species_subtree_time_seq = species_subtree.time_sequences(coalescent_process=species_subtree_coal_process)
        print('\nspecies_subtree_time_seq:')
        pprint.pprint(species_subtree_time_seq)

        # save subtree
        species_subtree.save_to_file(path='output/subtrees/species_subtree_' + str(node_id) + '_' + str(event['distance']*1000000)[:4])
        
        print('\ngene_subtree nodes:')
        gene_subtree = GeneTree(time_sequences=species_subtree_time_seq, species_tree=species_subtree)
        gene_subtree.print_nodes()

        print('\ngene_subtree dl_process:')
        gene_subtree_dup_events = gene_subtree.dl_process(lambda_dup=0.1, lambda_loss=0.03)
        print('\ngene_subtree dup_events:')
        print(gene_subtree_dup_events)

        self.duplication_subtree(coalescent_process=species_subtree_coal_process, dup_events=gene_subtree_dup_events)
        return

    def duplication_subtree(self, coalescent_process, dup_events):
        if (not dup_events):
            return
        for event in dup_events:
            node_id = None
            coal_distance = None
            for k, v in coalescent_process.items():
                for elem in v:
                    if (event['name'] in elem['to_set'] and event['name'] not in elem['from_set']):
                        node_id = int(k)
                        coal_distance = elem['distance']
            if (not node_id):
                print('to be considered...')
            else:
                self.__duplication_subtree_recurse(event=event, node_id=node_id, coal_distance=coal_distance)
        return