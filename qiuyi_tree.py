from io import StringIO
import skbio
import skbio.tree

import numpy as np
import random
import collections
import pprint

np.random.seed(6)


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
        self.nodes_dict = {}
        self.root = None
        self.leaves = []
        return

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
        # open table file and construct species tree
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
            self.nodes_dict[node.name] = node

        # find children
        for node in self.nodes:
            children, node.distance_to_children = self.children_distances(self.skbio_tree, node.name)
            for child in children:
                node_id = self.nodes_dict[child.name].node_id
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
                 newick_path,
                 lambda0):
        GenericTree.__init__(self)
        self.__construct_species_nodes(newick_path)
        for node in self.nodes:
            if (node.parent < 0):
                self.root = node
        self.leaves = [node.node_id for node in self.nodes if not node.children]

        self.set = [[str(node.node_id) + '*'] if not node.children else [] for node in self.nodes]
        self.lambda0 = lambda0
        self.coalescent_process = collections.defaultdict(list)
        return

    def __construct_species_nodes(self, newick_path):
        output_path = 'data/species_nodes_table.txt'
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
                             distance):
        if (len(self.set[node_id]) <= 1):
            return
        else:
            lambda_c = len(self.set[node_id]) * self.lambda0
            distance_fake = np.random.exponential(scale=1.0/lambda_c)
            if (distance < distance_fake):
                return
            else:
                if (len(self.set[node_id]) >= 2):
                    temp_set = sorted(self.set[node_id])
                    couple = np.random.choice(self.set[node_id], size=2, replace=False)
                    self.set[node_id] = [''.join(self.__star_sorted(couple))] + [e for e in self.set[node_id] if e not in couple]

                    # print process
                    print("initial node " + str(node_id) + ": " + str(temp_set))
                    print("coalescent at node " + str(node_id) + ": " + str(self.set[node_id]) + ", " + "distance = " + str(distance_fake))

                    # save process
                    self.coalescent_process[str(node_id)].append({
                        'from_set': temp_set, 
                        'to_set': self.set[node_id].copy(),
                        'distance': distance_fake
                    })
                else:
                    return
                distance = distance - distance_fake
                self.__coalescent_recurse(node_id=node_id, distance=distance)
        return

    def coalescent(self):
        old_leaves = self.leaves.copy()
        new_leaves = []
        labelled = [False for _ in range(len(self.nodes))]
        while (True):
            for leaf in old_leaves:
                if (leaf == self.root.node_id):
                    self.__coalescent_recurse(node_id=self.root.node_id, distance=10000)
                    break
                else:
                    parent = self.nodes[leaf].parent
                    children = self.nodes[parent].children
                    if (labelled[leaf]):
                        continue
                    labelled[leaf] = True
                    if (len(self.set[children[0]]) != 0 
                        and len(self.set[children[1]]) != 0):
                        self.__coalescent_recurse(node_id=children[0], 
                                                  distance=self.nodes[children[0]].distance_to_parent)
                        self.__coalescent_recurse(node_id=children[1], 
                                                  distance=self.nodes[children[1]].distance_to_parent)
                        self.set[parent] = list(set().union(self.set[children[0]], self.set[children[1]]))
                        if (len(new_leaves) > 0):
                            new_leaves = [e for e in new_leaves if e != children[0] and e != children[1]]
                        new_leaves.append(parent)
                    else:
                        new_leaves.append(leaf)
            if (leaf == self.root.node_id):
                break
            temp_new_leaves = []
            for new_leaf in new_leaves:
                if (new_leaf not in temp_new_leaves):
                    temp_new_leaves.append(new_leaf)
            old_leaves = temp_new_leaves.copy()
            new_leaves = []
            labelled = [False for _ in range(len(self.nodes))]

        print('\ncoalescent_process:')
        pprint.pprint(self.coalescent_process)
        return


class GeneTree(GenericTree):
    def __init__(self,
                 species_tree,
                 lambda_dup,
                 lambda_loss):
        GenericTree.__init__(self)
        self.leaves = species_tree.leaves
        
        self.species_nodes = species_tree.nodes
        self.coalescent_process = species_tree.coalescent_process
        self.total_distance = self.__distance_to_root_recurse(node_id=0)
        self.time_sequence = {}
        for leaf in self.leaves:
            self.time_sequence[str(leaf)] = self.__reverse_time_sequence(target_star=str(leaf)+'*')

        print('\ntime_sequence:')
        pprint.pprint(self.time_sequence)

        self.__construct_gene_nodes()
        for node in self.nodes:
            if (node.parent < 0):
                self.root = node
        
        self.lambda_dup = lambda_dup
        self.lambda_loss = lambda_loss

        return

    def __distance_to_root_recurse(self, node_id):
        if (self.species_nodes[node_id].parent < 0):
            return 0
        else:
            d2p = self.species_nodes[node_id].distance_to_parent
            parent = self.species_nodes[node_id].parent
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
    
    def __reverse_time_sequence(self, target_star):
        sequence = []
        for k, v in self.coalescent_process.items():
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
                            sequence += self.__reverse_time_sequence(target_star=e)
        return sequence

    def __star_replace(self, string, substring):
        a = string.split('*')[:-1]
        b = substring.split('*')[:-1]
        diff = set(a).difference(set(b))
        return ''.join([e + '*' for e in sorted(list(diff))])

    def __distance(self, node_name, parent_name):
        for leaf, sequence in self.time_sequence.items():
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

    def __construct_skbio_tree_recurse(self, skbio_tree_node):
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
        for leaf, sequence in self.time_sequence.items():
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
                    self.__construct_skbio_tree_recurse(child_one)
                    self.__construct_skbio_tree_recurse(child_two)
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
        tree.name = self.time_sequence['0'][-1][0]
        self.__construct_skbio_tree_recurse(tree)
        tree.length = None
        self.skbio_tree = tree
        super().newick_to_table(skbio_tree=tree, output_path='data/gene_nodes_table.txt')
        super().construct_nodes('data/gene_nodes_table.txt')

        return

    def __dl_process_recurse(self, tree, distance):
        node = self.nodes_dict[tree.name]
        distance_dup = np.random.exponential(scale=1.0/self.lambda_dup)
        distance_loss = np.random.exponential(scale=1.0/self.lambda_loss)
        if (distance_dup < distance_loss and distance_dup < distance):
            print('duplication at node ' + str(node.node_id) + ' (' + node.name + ')' + ' with distance ' + str(distance - distance_dup))
            self.__dl_process_recurse(tree, distance - distance_dup)
        elif (distance_loss <= distance_dup and distance_loss < distance):
            print('loss at node ' + str(node.node_id) + ' (' + node.name + ')' + ' with distance ' + str(distance_loss))
        else:
            print('nothing happened at node ' + str(node.node_id) + ' (' + node.name + ')')
            if (node.children):
                child_one = tree.children[0]
                child_two = tree.children[1]
                distance_to_child_one = node.distance_to_children[0]
                distance_to_child_two = node.distance_to_children[1]
                self.__dl_process_recurse(child_one, distance_to_child_one)
                self.__dl_process_recurse(child_two, distance_to_child_two)
            else:
                print('reach the end of node ' + str(node.node_id) + ' (' + node.name + ')')

    def dl_process(self):
        self.__dl_process_recurse(self.skbio_tree, distance=0)
        
