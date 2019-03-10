from io import StringIO
from skbio import read
from skbio.tree import TreeNode

import qiuyi_tree


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

def main():
    qtree = qiuyi_tree.SpeciesTree(table_file_path='data/nodes_table.txt',
                                   lambda0=1)
    qtree.print_nodes()
    qtree.coalescent()


if __name__ == "__main__":
    main()
