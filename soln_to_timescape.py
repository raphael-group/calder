# Author: Matt Myers (created 6/12/2019)

from string import ascii_uppercase
import pandas as pd
import networkx as nx
import sys
import os
import traceback
from math import log

def to_ascii(x):
    """ Converts a 0-indexed integer to the corresponding Excel-style column name """
    x += 1
    result = []
    while x:
        x, rem = divmod(x-1, 26)
        result[:0] = ascii_uppercase[rem]
    return ''.join(result)

def parse_csv(csv_fname):
    # Read in solution CSV file (names=range(10000) specifies the max number of mutations))
    df = pd.read_csv(csv_fname, header=None, names=range(10000))

    # Separate rows of CSV file into 2 tables
    table_names = ["Fhat", "U"]
    groups = df[0].isin(table_names).cumsum()
    tables = {g.iloc[0,0]: g.iloc[1:] for k,g in df.groupby(groups)}

    # Drop empty columns
    tables['Fhat'].dropna(axis=1, inplace=True)
    tables['U'].dropna(axis=1, inplace=True)

    # Set up header rows
    for key, table in tables.items():
        header = table.iloc[0]
        df = table[1:]
        df.columns = header
        tables[key] = df
    return tables

def parse_dot(dot_fname):
    # Read in DOT file
    tree = nx.drawing.nx_pydot.read_dot(dot_fname)

    # Extract vertex labels
    vertex_labels = {k:tree.nodes[k]['label'].strip("\"") for k in tree.nodes}

    # Initialize tree dictionary using vertex labels
    # {source1: [destination1, destination2, ...], source2: [...], ...}
    my_edges = [(vertex_labels[x], vertex_labels[y]) for (x, y, _) in tree.edges]
    tree_dict = {}
    for (a, b) in my_edges:
        if a in tree_dict:
            tree_dict[a].append(b)
        else:
            tree_dict[a] = [b]

    return tree_dict

def main():
    try:
        tables = parse_csv(sys.argv[1])
    except Exception as e:
        print("Error parsing argument 1 as CSV file.")
        traceback.print_exc()
        quit(1)
    
    try:
        tree_dict = parse_dot(sys.argv[2])
    except Exception as e:
        print("Error parsing argument 2 as DOT file.")
        traceback.print_exc()
        quit(1)

    try:
        write_timescape(tables['U'], tree_dict, sys.argv[3])
    except Exception as e:
        print("Error converting file to Timescape format.")
        traceback.print_exc()
        quit(1)
    

def write_timescape(U, tree_dict, stem):
    rows = []
    rows.append(["timepoint", "clone_id", "clonal_prev"])
    
    # (topological sorting is not necessary, but makes labels more intuitive)
    order = topological_sort(tree_dict)
    indexes = {}
    for i in range(len(order)):
        indexes[order[i]] = i    

    m, n = U.shape
    n -= 1 # first column contains sample labels
    for t in range(m):
        tp = "T%d" % (t + 1)
        for i in range(n):
            newrow = []
            newrow.append(tp)
            newrow.append(to_ascii(indexes[U.columns[i + 1]]))
            newrow.append(str(U.iloc[t][i + 1]))
            rows.append(newrow)
    
    with open(stem + "_prev.txt", "w") as f:
        for outrow in rows[:-1]:
            f.write("\t".join(outrow) + os.linesep)
        f.write('\t'.join(rows[-1]))
    
    with open(stem + "_edges.txt", "w") as f:
        f.write("\t".join(['source', 'target']))
        f.write(os.linesep)
        for source, dests in tree_dict.items():
            for dest in dests:
                f.write("\t".join([to_ascii(indexes[source]), to_ascii(indexes[dest])]))
                f.write(os.linesep)

    print("Vertex labels are abbreviated in Timescape files using the following mapping:")
    for i in range(len(order)):
        print("Clone %s: %s" % (to_ascii(i), order[i]))

def topological_sort(tree):
    roots = {}
    for a in tree.keys():
        roots[a] = True
    for _, dests in tree.items():
        for dest in dests:
            if dest in roots:
                del roots[dest]
    assert len(roots) == 1
    root = list(roots.keys())[0]
    
    S = [root]
    L = []
    while len(S) > 0:
        curr = S[0]
        S = S[1:]
        L.append(curr)
        if curr in tree:
            for dest in tree[curr]:
                S.append(dest)
    return L

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python soln_to_timescape.py [CSV solution file] [DOT solution file] [output filename stem]")
        quit(1)
    else:
        main()