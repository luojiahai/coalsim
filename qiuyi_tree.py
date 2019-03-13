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
                 node_id,
                 parent,
                 distance_to_parent):
        self.node_id = node_id
        self.parent = parent
        self.distance_to_parent = distance_to_parent
        self.childs = []
        self.distance_to_childs = []
        return
    
    def __repr__(self):
        return "node_id: {}, parent: {}, distance_to_parent: {}, childs: {}".format(self.node_id, 
                                                                        self.parent, 
                                                                        self.distance_to_parent,
                                                                        self.childs)


class SpeciesTree(object):
    def __init__(self,
                 table_file_path,
                 lambda0):
        self.lambda0 = lambda0

        self.nodes = []
        self.__construct_species_nodes(table_file_path)

        self.set = [[str(node.node_id) + '*'] if not node.childs else [] for node in self.nodes]
        self.leaves = [node.node_id for node in self.nodes if not node.childs]
        self.root = None
        for node in self.nodes:
            if (node.parent < 0):
                self.root = node

        self.coalescent_process = collections.defaultdict(list)
        return

    def print_nodes(self):
        for node in self.nodes:
            print(node)
        return

    def __construct_species_nodes(self, path):
        # open table file and construct species tree
        f = open(path, 'r')
        f.readline()
        for line in f:
            splited = line.strip().split('\t')
            parent = splited[2]
            if (parent == 'None'):
                parent = -1
            d2p = splited[3]
            if (d2p == 'None'):
                d2p = -1.0
            self.nodes.append(TreeNode(node_id=int(splited[0]),
                                       parent=int(parent),
                                       distance_to_parent=float(d2p)))
        
        # find childs
        childs = [[] for _ in range(len(self.nodes))]
        for node in self.nodes:
            if (node.parent < 0): continue
            childs[node.parent].append(node.node_id)
        for node in self.nodes:
            node.childs = sorted(childs[node.node_id])

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
                    childs = self.nodes[parent].childs
                    if (labelled[leaf]):
                        continue
                    labelled[leaf] = True
                    if (len(self.set[childs[0]]) != 0 
                        and len(self.set[childs[1]]) != 0):
                        self.__coalescent_recurse(node_id=childs[0], distance=self.nodes[childs[0]].distance_to_parent)
                        self.__coalescent_recurse(node_id=childs[1], distance=self.nodes[childs[1]].distance_to_parent)
                        self.set[parent] = list(set().union(self.set[childs[0]], self.set[childs[1]]))
                        if (len(new_leaves) > 0):
                            new_leaves = [e for e in new_leaves if e != childs[0] and e != childs[1]]
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

        pprint.pprint(self.coalescent_process)
        return


class GeneTree(object):
    def __init__(self,
                 species_tree):
        self.species_nodes = species_tree.nodes
        self.coalescent_process = species_tree.coalescent_process
        self.leaves = species_tree.leaves
        self.total_distance = self.__distance_to_root_recurse(node_id=0)

        self.time_sequence = {}
        for leaf in self.leaves:
            self.time_sequence[leaf] = self.__reverse_time_sequence(target_star=str(leaf)+'*')
        pprint.pprint(self.time_sequence)

        self.nodes = []
        self.__construct_gene_nodes()
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

    def __in_set(self, target, clade):
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
                        if self.__in_set(target_star, e):
                            couple = e.replace(target_star, '')
                            walking_distance = self.__walking_distance(int(k), branch_distance=branch_distance)
                            # pair = (couple, walking_distance)
                            pair = (e, walking_distance)
                            sequence.append(pair)
                            sequence += self.__reverse_time_sequence(target_star=e)
        return sequence

    def __construct_gene_nodes(self):
        return


def newick_to_table(output_path, input_path=None, skbio_tree=None):
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
            ret = sorted(
                ret, key=lambda x: x['distance_to_root'], reverse=True)
        return ret

    def _output_to_file(path, nodes):
        f = open(path, 'w')
        f.write('id\tname\tparent\td2p\n')
        for i in range(len(nodes)):
            parent_index = 'None'
            for j in range(len(nodes)):
                if (nodes[i]['parent'] is nodes[j]['object']):
                    parent_index = str(j)
            f.write(str(i) + '\t' + nodes[i]['name'] + '\t' +
                    parent_index + '\t' + str(nodes[i]['distance']) + '\n')
        f.close()

    root = _parse(tree)
    root = _rename(root)
    nodes = _to_list(root, root=root['object'])
    _output_to_file(output_path, nodes)