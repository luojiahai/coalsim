from io import StringIO
from skbio import read
from skbio.tree import TreeNode
import numpy as np
import random


def newick_to_table():
    f = open('data/tree_sample.txt')
    tree = read(f, format="newick", into=TreeNode)
    f.close()

    def parse(tree):
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
            node['children'].append(parse(children))
        return node

    def rename(node):
        ret = node.copy()
        name = ''
        for i in range(len(ret['children'])):
            ret['children'][i] = rename(ret['children'][i])
            name += ret['children'][i]['name']
        if (not ret['name']):
            ret['name'] = name
        return ret

    def to_list(node, root):
        d = node.copy()
        ret = []
        for i in range(len(d['children'])):
            ret += to_list(d['children'][i], root=root)
        del d['children']
        d['distance_to_root'] = d['object'].distance(root)
        ret.append(d)
        if (d['object'] is root):
            ret = sorted(
                ret, key=lambda x: x['distance_to_root'], reverse=True)
        return ret

    def output_to_file(path, nodes):
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

    root = parse(tree)
    root = rename(root)
    nodes = to_list(root, root=root['object'])
    output_to_file('data/nodes_table.txt', nodes)


def read_table(path):
    f = open(path, 'r')
    f.readline()
    temp = []
    parents = []
    distances = []
    names = []
    for line in f:
        splited = line.strip().split('\t')
        parent = splited[2]
        if (parent == 'None'):
            parent = -1
        d2p = splited[3]
        if (d2p == 'None'):
            d2p = 0.0
        temp.append({
            'id': int(splited[0]),
            'parent': int(parent)
        })
        parents.append(int(parent))
        distances.append(float(d2p))
        names.append(splited[1])
    childs = [[] for _ in range(len(temp))]
    for e in temp:
        if (e['parent'] < 0):
            continue
        childs[e['parent']].append(e['id'])
    return 

def coal(node, mark, distance, set1, lambda0):
    if (mark <= 1):
        return (mark, set1, distance)
    else:
        lambda_c = mark * lambda0
        distance_fake = np.random.exponential(scale=1.0/lambda_c)
        if (distance < distance_fake):
            return (mark, set1, distance)
        else:
            # change set
            if (len(set1[node]) >= 2):
                print("initial node " + str(node) + ": " + str(set1[node]))
                couple = random.sample(set1[node], 2)
                set1[node] = [int(''.join([str(e) for e in sorted(couple)]))] + [e for e in set1[node] if e not in couple]
                print("coalescent at node " + str(node) + ": " + str(set1[node]) + ", " + "distance = " + distance_fake)
            else:
                return (mark, set1, distance)
            distance = distance - distance_fake
            coal(node, mark - 1, distance, set1, lambda0)

def main():
    # newick_to_table()
    data = read_table('data/nodes_table.txt')

    print(data)
    # n_nodes = data['names']
    # mark = np.zeros(10)
    


if __name__ == "__main__":
    main()
