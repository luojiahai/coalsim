import numpy as np
import random


class TreeNode(object):
    def __init__(self,
                 node_id,
                 parent,
                 distance_to_parent):
        self.node_id = node_id
        self.parent = parent
        self.distance_to_parent = distance_to_parent
        self.childs = []
        self.walking_distance = -1
    
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
        self.set = [[str(node.node_id)] if not node.childs else [] for node in self.nodes]
        self.leaves = [node.node_id for node in self.nodes if not node.childs]
        self.root = None
        for node in self.nodes:
            if (node.parent < 0):
                self.root = node

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
                    print("initial node " + str(node_id) + ": " + str(self.set[node_id]))
                    couple = random.sample(self.set[node_id], 2)
                    self.set[node_id] = [''.join([str(e) for e in couple])] + [str(e) for e in self.set[node_id] if e not in couple]
                    print("coalescent at node " + str(node_id) + ": " + str(self.set[node_id]) + ", " + "distance = " + str(distance_fake))
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
            new_leaves = temp_new_leaves
            self.leaves = new_leaves.copy()
            labelled = [False for _ in range(len(self.nodes))]
            new_leaves = []
