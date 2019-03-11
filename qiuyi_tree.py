import numpy as np
import random
import collections

np.random.seed(3)


class TreeNode(object):
    def __init__(self,
                 node_id,
                 parent,
                 distance_to_parent):
        self.node_id = node_id
        self.parent = parent
        self.distance_to_parent = distance_to_parent
        self.childs = []
    
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
        self.__table_to_nodes(table_file_path)
        self.set = [[str(node.node_id) + '*'] if not node.childs else [] for node in self.nodes]
        self.leaves = [node.node_id for node in self.nodes if not node.childs]
        self.root = None
        for node in self.nodes:
            if (node.parent < 0):
                self.root = node
        self.coalescent_process = collections.defaultdict(list)

    def __table_to_nodes(self, path):
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
        self.__find_childs()

    def __find_childs(self):
        childs = [[] for _ in range(len(self.nodes))]
        for node in self.nodes:
            if (node.parent < 0): continue
            childs[node.parent].append(node.node_id)
        for node in self.nodes:
            node.childs = sorted(childs[node.node_id])
                                    
    def print_nodes(self):
        for node in self.nodes:
            print(node)

    def __star_sorted(self, couple):
        string = ''
        for e in couple:
            string += e
        splited = string.split('*')[:-1]
        splited = sorted([int(e) for e in splited])
        return [str(e) + '*' for e in splited]

    def __coal_recurse(self,
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
                    print("initial node " + str(node_id) + ": " + str(temp_set))
                    couple = random.sample(self.set[node_id], 2)
                    self.set[node_id] = [''.join(self.__star_sorted(couple))] + [e for e in self.set[node_id] if e not in couple]
                    print("coalescent at node " + str(node_id) + ": " + str(self.set[node_id]) + ", " + "distance = " + str(distance_fake))
                    self.coalescent_process[str(node_id)].append({
                        'from_set': temp_set, 
                        'to_set': self.set[node_id].copy(),
                        'distance': distance_fake
                    })
                else:
                    return
                distance = distance - distance_fake
                self.__coal_recurse(node_id=node_id, distance=distance)

    def coalescent(self):
        new_leaves = []
        labelled = [False for _ in range(len(self.nodes))]
        while (True):
            for leaf in self.leaves:
                if (leaf == self.root.node_id):
                    self.__coal_recurse(node_id=self.root.node_id, distance=10000)
                    break
                else:
                    parent = self.nodes[leaf].parent
                    childs = self.nodes[parent].childs
                    if (labelled[leaf]):
                        continue
                    labelled[leaf] = True
                    if (len(self.set[childs[0]]) != 0 
                        and len(self.set[childs[1]]) != 0):
                        self.__coal_recurse(node_id=childs[0], distance=self.nodes[childs[0]].distance_to_parent)
                        self.__coal_recurse(node_id=childs[1], distance=self.nodes[childs[1]].distance_to_parent)
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
            self.leaves = temp_new_leaves.copy()

            new_leaves = []
            labelled = [False for _ in range(len(self.nodes))]


class GeneTree(object):
    def __init__(self,
                 species_tree):
        self.nodes = species_tree.nodes
        self.coalescent_process = species_tree.coalescent_process
        self.leaves = [node.node_id for node in self.nodes if not node.childs]
        self.total_distance = self.__distance_to_root(node_id=0)

        for leaf in self.leaves:
            print(leaf)
            self.__extract_couples(target_star=str(leaf)+'*')

    def __distance_to_root(self, node_id):
        if (self.nodes[node_id].parent < 0):
            return 0
        else:
            d2p = self.nodes[node_id].distance_to_parent
            parent = self.nodes[node_id].parent
            return d2p + self.__distance_to_root(parent)

    def __in_set(self, target, clade):
        if (len(target) < len(clade)):
            splited_target = target.split('*')[:-1]
            splited_clade = clade.split('*')[:-1]
            return set(splited_target).issubset(set(splited_clade))
        else:
            return False
    
    def __extract_couples(self, target_star):
        sequence = []
        for k, v in self.coalescent_process.items():
            for elem in v:
                if (target_star in elem['from_set'] and target_star not in elem['to_set']):
                    for e in elem['to_set']:
                        if self.__in_set(target_star, e):
                            couple = e.replace(target_star, '')
                            walking_distance = elem['distance'] + (self.total_distance - self.__distance_to_root(int(k)))
                            # pair = (couple, walking_distance)
                            pair = (e, walking_distance)
                            print(pair)
                            self.__extract_couples(target_star=e)
        return sequence


