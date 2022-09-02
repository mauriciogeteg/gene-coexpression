import os, sys, time, random
import networkx as nx, pandas as pd, numpy as np
from collections import Counter
import pickle as pk

def read_expression(file, sep=None):
    names = []
    mat = []
    with open(file, "r") as f:
        headers = f.readline()
        for line in f.readlines():
            name, tmp = line.strip().split(sep, 1)
            names.append(name)
            mat.append(np.array([float(x) for x in tmp.strip().split(sep)]))
    return names, mat, headers

if __name__ == '__main__':
    print(time.strftime("%a, %d %b %Y %H:%M:%S"))
    file_expression = "expression_data/sugar_cane_expression_2018.csv"
    file_mapping1 = "blastn_2019_2018"
    file_mapping2 = "blastn_2018_2019"
    file_result = "expression_data/sugar_cane_allele_2019.csv"
    df_map = pd.read_csv(file_mapping1, sep=None, index_col=None, names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"], engine='python')
    names, exp, headers = read_expression(file_expression, ",")
    #print(df_map.head())
    #print(df_map[df_map['pident']<100.0])
    #print(df_map[df_map.duplicated('sseqid')])
    # take duplicated queries
    df_dup = df_map[df_map.duplicated('sseqid')]
    rep_CDS = df_dup['sseqid'].unique()
    tmp = df_map[df_map['sseqid'].isin(rep_CDS)].sort_values(by=['sseqid'])
    #print(tmp)

    # Bipartite graph creation
    B = nx.Graph()
    B.add_nodes_from(df_map['qseqid'].tolist(), bipartite=0)
    B.add_nodes_from(df_map['sseqid'].tolist(), bipartite=1)
    for index,row in df_map.iterrows():
        if not B.has_edge(row['qseqid'], row['sseqid']):
            B.add_edge(row['qseqid'], row['sseqid'], capacity=1, weight=int(row['pident']*-1000))
        elif B[row['qseqid']][row['sseqid']]['weight']>int(row['pident']*-1000):
            B[row['qseqid']][row['sseqid']]['weight'] = int(row['pident']*-1000)
    print("2019->2018", row['qseqid'], row['sseqid'])
    print("Nodes", B.number_of_nodes(), ", Edges", B.number_of_edges())

    # now load the opposite direction, but change the first two column names
    df_map = pd.read_csv(file_mapping2, sep=None, index_col=None, names=["sseqid", "qseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"], engine='python')
    B.add_nodes_from(df_map['qseqid'].tolist(), bipartite=0)
    B.add_nodes_from(df_map['sseqid'].tolist(), bipartite=1)
    for index,row in df_map.iterrows():
        if not B.has_edge(row['qseqid'], row['sseqid']):
            B.add_edge(row['qseqid'], row['sseqid'], capacity=1, weight=int(row['pident']*-1000))
        elif B[row['qseqid']][row['sseqid']]['weight']>int(row['pident']*-1000):
            B[row['qseqid']][row['sseqid']]['weight'] = int(row['pident']*-1000)
    print("2018->2019", row['qseqid'], row['sseqid'])
    print("Nodes", B.number_of_nodes(), ", Edges", B.number_of_edges())
    df_map = pd.read_csv(file_mapping1, sep=None, index_col=None, names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"], engine='python')
    
    # maximum matching of each component of B
    # (computed only when the number of nodes to the right is >1)
    restricted = dict()
    ct = 0
    l_nodes = set(df_map['qseqid'].tolist())
    r_nodes = set(B) - l_nodes
    # CONFIGURING SOURCE AND TARGET
    B.add_node('SS')
    B.add_node('ST')
    for node in r_nodes:
        B.add_edge('SS', node, capacity=1)  # zero weight
    for node in l_nodes:
        B.add_edge(node, 'ST')  # infinite capacity, zero weight
    maxflow = nx.max_flow_min_cost(B, 'SS', 'ST')
    for v in B['SS']:
        if maxflow['SS'][v] == 0: continue
        ct += 1
        for u in B[v]:
            if u == 'SS': continue
            if maxflow[v][u] > 0 and maxflow[u][v]==0:
                if v in restricted:
                    raise Exception(f"{v} appears assigned twice! {u} and {restricted[v]}")
                else:
                    restricted[v] = u
    print("Number of transcripts:", ct)
    
    # # computing maximal bipartite matching
    # for comp in nx.connected_components(B):
        # sub = B.subgraph(comp)
        # if sub.number_of_nodes() > 2: ct += 1
        # test = Counter([x in r_nodes for x in sub.nodes()])
        # if test[True] == 1: continue
        # matches = nx.bipartite.maximum_matching(sub)
        # #print(sub.edges())
        # #print("matches", matches)
        # for key in matches:
            # if key in r_nodes:
                # restricted[key] = matches[key]
    # print("Total of non-trivial components:", ct)
    # print("Number of restricted connections:", len(restricted))

    new_names = sorted(df_map['qseqid'].unique().tolist())
    dict_old = {v:u for u,v in enumerate(names)}
    dict_new = {v:u for u,v in enumerate(new_names)}
    mat = np.zeros( (len(new_names), len(exp[0])) )
    print("Number of genes:", len(new_names))
    for key in restricted:
        value = restricted[key]
        mat[dict_new[value]] += exp[dict_old[key]]
        
    # for index,row in df_map.iterrows():
        # if row['sseqid'] not in restricted:
            # #print(index)
            # #print(row['sseqid'])
            # #print(dict_old[row['sseqid']])
            # mat[dict_new[row['qseqid']]] += exp[dict_old[row['sseqid']]]
        # elif row['qseqid'] == restricted[row['sseqid']]:
            # mat[dict_new[row['qseqid']]] += exp[dict_old[row['sseqid']]]

    # Saving the resulting DataFrame
    with open(file_result, "w") as f:
        f.write(headers)
        for i,name in enumerate(new_names):
            f.write(name)
            for x in mat[i]:
                f.write("," + str(float(x)))
            f.write("\n")
    pk.dump(restricted, open("restricted.pk", "wb"))
    with open("restricted.txt", "w") as f:
        f.write("key,value\n")
        for key in restricted:
            f.write(f"{key},{restricted[key]}\n")
    
    print("DONE!")
