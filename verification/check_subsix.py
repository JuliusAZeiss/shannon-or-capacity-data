#!/usr/bin/env python3
"""Exact finite checks for the compact sub-six proof; standard library only.

The matrix data are in compact_bd.py and compact_h.py. No solver or floating
point arithmetic is used. The independent source-derived checker separately
checks the manuscript tables. This program certifies the finite inputs, not
the analytical spectral interpolation theorem.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import hashlib, importlib.util, json

def load(name):
    path=Path(__file__).resolve().parent/(name+'.py')
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

bd=load('compact_bd')
hc=load('compact_h')

H_ROWS = ('010101110001100100011101010', '100111111101001000001100000', '000000011010011000101100101', '110011110100101010101101001', '010101101011100111011010000', '110110011100001010100001110', '110110000110110011100011111', '111101001101000100000010000', '011011010111101100101101100', '010101111010100000001010001', '001010101101111000001111000', '110010011010111110110010011', '100110101111001100110010000', '001000100011001110100101111', '011101001011110010001000101', '100010011001110011101011000', '000111100001011101101010100', '000010100000000110111101011', '001101101001110111010110110', '100010000001100001100101110', '111110001110001111000011011', '111100001010010001110011010', '000010110111100110101101010', '100101101010010101011110011', '001001101000011010110000001', '100001100001010001111111000', '001100100101011001001001100')

COLOR_CLASSES = ((1, 2, 12, 16, 23), (2, 5, 6, 15, 19), (7, 10, 16, 25, 26), (7, 12, 20, 21, 24), (0, 8, 13, 17, 22), (3, 8, 13, 19, 22), (1, 10, 15, 25, 26), (0, 4, 9, 14, 18), (4, 9, 14, 18, 23), (5, 6, 11, 20, 21), (3, 11, 17, 24))

def require(condition, message):
    if not condition:
        raise AssertionError(message)

def dot3(u, v):
    return sum(x * y for x, y in zip(u, v)) % 3

def rank3(a):
    """Independent Gaussian elimination over F3, for the fitting controls."""
    m = [[x % 3 for x in row] for row in a]
    nrows = len(m)
    ncols = len(m[0]) if m else 0
    require(all(len(row) == ncols for row in m), "ragged matrix")
    rank = 0
    for col in range(ncols):
        pivot = next((i for i in range(rank, nrows) if m[i][col]), None)
        if pivot is None:
            continue
        m[rank], m[pivot] = m[pivot], m[rank]
        inv = pow(m[rank][col], -1, 3)
        m[rank] = [(inv * x) % 3 for x in m[rank]]
        for i in range(rank + 1, nrows):
            scale = m[i][col]
            if scale:
                m[i] = [(x - scale * y) % 3 for x, y in zip(m[i], m[rank])]
        rank += 1
        if rank == nrows:
            break
    return rank

def square_code_check(a, code, size):
    require(len(code) == size and len(set(code)) == size, "wrong code size or repeated word")
    require(all(len(w) == 2 and all(0 <= x < len(a) for x in w) for w in code), "invalid code letter")
    require(all(a[u[0]][v[0]] or a[u[1]][v[1]] for u, v in combinations(code, 2)),
            "square code has a nonadjacent pair")
    return size * (size - 1) // 2

def graph_check(a):
    n = len(a)
    require(all(len(row) == n for row in a), "nonsquare adjacency")
    require(all(a[i][i] == 0 for i in range(n)), "graph has a loop")
    require(all(a[i][j] in (0, 1) and a[i][j] == a[j][i]
                for i in range(n) for j in range(n)), "graph is not simple undirected")

def graph_digest(a):
    bits = "\n".join("".join(str(int(x)) for x in row) for row in a)
    return hashlib.sha256(bits.encode("ascii")).hexdigest()

def main():
    vectors=bd.VECTORS
    require(len(vectors)==27 and len(set(vectors))==27,"source order")
    require(all(dot3(v,v) for v in vectors),"isotropic source vector")
    B,D=bd.adjacency(vectors)
    H=[[int(x) for x in row] for row in H_ROWS]
    G=[[0]*54 for _ in range(54)]
    for i in range(27):
        for j in range(27):G[i][j]=B[i][j];G[i+27][j+27]=D[i][j]
    for graph in B,D,H,G:graph_check(graph)
    scalar=[[pow(dot3(u,u),-1,3)*dot3(u,w)%3 for w in vectors] for u in vectors]
    require(all(scalar[i][j]==int(i==j) for i in range(27) for j in range(27)
                if i==j or B[i][j]),"scalar fitting condition")
    require(rank3(scalar)==5,"rank five")
    labels=[v for v in range(29) if v not in (1,17)];pos={v:i for i,v in enumerate(labels)}
    require(len(H)==27 and sum(map(sum,H))==352,"27-vertex, 176-edge partner")
    h_code=[(i,12*i%29) for i in range(29) if i not in (1,17,28)]+[(28,12)]
    counts=[square_code_check(G,[(i,i+27) for i in range(27)],27),
            square_code_check(H,[(pos[x],pos[y]) for x,y in h_code],27)]
    require(len(COLOR_CLASSES)==11,"eleven colours")
    for colour in COLOR_CLASSES:
        require(all(not H[i][j] for i,j in combinations(colour,2)),"independent class")
    assigned=[[a for a,colour in enumerate(COLOR_CLASSES) if i in colour] for i in range(len(labels))]
    require(all(len(x)==2 for x in assigned),"twofold colouring")
    basis=[[int(k==a) for k in range(11)] for pair in assigned for a in pair]
    fit=[[sum(a*b for a,b in zip(u,w)) for w in basis] for u in basis]
    require(rank3(fit)==11,"rank eleven")
    for i in range(len(labels)):
        for j in range(len(labels)):
            if i==j or H[i][j]:
                require(all(fit[2*i+a][2*j+b]==int(i==j and a==b)
                            for a in range(2) for b in range(2)),"block fitting condition")
    sources=bd.verify();partner=hc.verify_compact_h(H,labels)
    require(hc.L**2<27 and max(F(16,3),hc.B,hc.Q)<6,"capacity comparisons")
    digest=lambda matrix:hashlib.sha256(json.dumps(matrix,separators=(',',':')).encode()).hexdigest()
    return {'status':'PASS: compact exact 54-by-27 certificate',
            'graph_orders':{'G':54,'H':27},
            'graph_SHA256':{'B':graph_digest(B),'H':graph_digest(H)},
            'square_code_pairs_checked':counts,'F3_fitting_ranks':{'B_scalar':5,'H_block_two':11},
            'bounds':{'B':str(F(16,3)),'D':str(hc.B),'H':str(hc.Q)},
            'matrix_SHA256':{'B':sources['B']['matrix_sha256'],'D':sources['D']['matrix_sha256'],
                             'H':digest(partner['matrix'])},
            'free_pair_orbits':{'B':23,'D':19,'H':91},
            'source_group_order':sources['group_order'],
            'source_blocks':{k:sources[k]['blocks'] for k in ('B','D')},
            'partner_block_orders':[len(partner['S']),len(partner['T'])],
            'positive_full_matrix_pivots':{'B':27,'D':27,'H':27},
            'interpolation_weight':str(hc.WEIGHT),'strict_power_comparisons':[True,True],
            'strict_product_upper':str(hc.L),'threshold_interval':'[5.193957,sqrt(27))',
            'scope':'Exact finite inputs; analytical spectral argument is proved in the supplement.'}

if __name__=='__main__':
    print(json.dumps(main(),indent=2))
