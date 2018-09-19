import re
import pandas as pd

from ete3 import Tree, TreeStyle, NodeStyle, TextFace
from HapFileIO import tree_info
from DataFunc import hap_made

hap_data = pd.read_table('/Users/zhengfuni/Desktop/Research/philippine/Y/1000Genome/1KG.hap', sep='\t', header=0, encoding='utf-8')
hap_type, filter, classify = ('key', True, 'smart')

def build_base_tree():
    base = ''
    level,length = 0, 0.5
    trunk_list = []
    trunk_length = {}
    count_num = tree_info[-1].count('-')
    for i in tree_info[1:]:
        level_new = i.count('-')
        trunk = i.split('\t')[-1]
        trunk_list.append(trunk)
        dif = level_new - level
        if dif == 0:
            length = length
            base = trunk + ',' + base
            if tree_info.index(i)+1 == len(tree_info):
                base =  '('*(base.count(')')-base.count('(')) + base
        elif dif > 0:
            length = 1+1*(count_num-level_new)
            base = trunk + ')' + base
            if tree_info.index(i)+1 == len(tree_info):
                base =  '('*(base.count(')')-base.count('(')) + base
        elif dif < 0:
            length = 1+1*(count_num-level_new)
            if tree_info.index(i)+1 != len(tree_info):
                base =  trunk + ',' + '('*abs(dif) + base
            else:
                base =  '('*(base.count(')')-base.count('(')) + trunk + base
        trunk_length[trunk] = length
        level = level_new
    base = '(' + base + 'Y-Adam);'
    base_tree = Tree(base, format=1)

    for node in base_tree.iter_descendants('postorder'):
        trunk = base_tree&node.name
        trunk.add_face(TextFace(node.name), column=0, position='branch-top')
        if node.is_leaf():
            node.dist = trunk_length[node.name]
        else:
            node.dist = 0.95

    ns = NodeStyle()
    ns['shape'] = 'circle'
    ns['size'] = 2
    #ns['hz_line_width'] = 1
    #ns['vt_line_width'] = 1
    ns['fgcolor'] = 'black'
    for n in base_tree.traverse():
        n.set_style(ns)

    #base_tree.show(tree_style=ts)
    return base_tree, trunk_list

def process_info(hap_data, hap_type, filter, classify, base_tree, trunk_list):
    hap_df = hap_made(hap_data.iloc[:, 0:5], hap_type, filter)
    initial = set(hap_df['Haplogroup'].map(lambda x: x[0]))
    delete_trunk = [x for x in trunk_list if x not in initial]

    '''
    for i in delete_trunk:
        node = base_tree.search_nodes(name=i)[0]
        node.delete()
    '''

    ts = TreeStyle()
    ts.show_leaf_name = False
    #ts.show_branch_length = True
    #ts.scale = 50
    ts.branch_vertical_margin = 0
    ts.optimal_scale_level = 'full'
    ts.show_scale = False
    base_tree.render('/Users/zhengfuni/Desktop/mytree2.png', tree_style=ts)
    #print(base_tree.get_ascii(show_internal=True))
    '''
    for i in initial:
        initial_group = sorted(list(filterfalse(lambda x: x[0] != i, hap_df['Haplogroup'])), key=lambda x: len(x))
    if classify == 'rough':
        hap_df['Haplogroup'] = hap_df['Haplogroup'].map(lambda x: x[0])
        hap_df['Mutation'] = [key_info.at[key_info['haplogroup'].tolist().index(i), 'mutation'] for i in hap_df['Haplogroup']]
    elif classify == 'smart':
        initial = set(hap_df['Haplogroup'].map(lambda x: x[0]))
        common_list = []
        for i in initial:
            initial_group = sorted(list(filterfalse(lambda x: x[0] != i, hap_df['Haplogroup'])), key=lambda x: len(x))
            while True:
                most_common = Counter(initial_group).most_common(1)[0][0]
                initial_group = [most_common if x.startswith(most_common) and len(x) > len(most_common) else x for x in initial_group]
                hap_num = len(most_common)-1
                for j in range(len(most_common)-1):
                    hap = most_common[0:hap_num]
                    if hap in initial_group:
                        initial_group = [hap if x == most_common else x for x in initial_group]
                        most_common = hap
                    hap_num-=1
                common_list.append(most_common)
                initial_group = list(filterfalse(lambda x: x.startswith(most_common), initial_group))
                if len(initial_group) == 0:
                    break
    '''



def main():
    base_tree, trunk_list = build_base_tree()
    process_info(hap_data, hap_type, filter, classify, base_tree, trunk_list)

if __name__ == '__main__':
    main()
