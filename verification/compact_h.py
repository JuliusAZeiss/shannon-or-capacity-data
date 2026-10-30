#!/usr/bin/env python3
"""Exact compact Gram certificate for the manuscript's 27-vertex H.

Run `python3 compact_h.py` from any directory. The graph is read from the
adjacency table in ../subsix_supplement.tex. There is no solver dependency.
The 91 triples below give one nonedge entry per orbit of the involution
on H. Its 27-by-27 matrix is positive definite, certified both by its
15- and 12-dimensional symmetry blocks and by direct rational elimination.
"""
from fractions import Fraction
from itertools import combinations
from pathlib import Path
import json
import re

SCALE=100000000
DIAGONAL=419656120
Q=Fraction(12991403,2500000)
B=Fraction(32437757,6250000)
L=Fraction(5193957,1000000)
WEIGHT=Fraction(13,1000)
H_LABELS=tuple(v for v in range(29) if v not in (1,17))
H_ORBIT_ENTRIES = (
    (0, 3, -92182906),
    (0, 5, 112385240),
    (0, 9, 105366413),
    (0, 10, 260048876),
    (0, 11, -64530825),
    (0, 14, 103723987),
    (0, 28, 24358075),
    (2, 3, 66703246),
    (2, 11, 207795723),
    (2, 13, 104538630),
    (2, 14, -136454182),
    (2, 16, 177512886),
    (2, 18, 115342751),
    (2, 19, -34383691),
    (2, 20, -100268450),
    (2, 21, 76229467),
    (2, 24, -95624150),
    (2, 25, 120690256),
    (2, 26, 152316144),
    (2, 27, 18986268),
    (2, 28, 103951980),
    (3, 4, -35794660),
    (3, 5, -98567795),
    (3, 6, 149505933),
    (3, 7, 190124261),
    (3, 10, -22994873),
    (3, 12, 66710199),
    (3, 13, 124899718),
    (3, 16, 61341650),
    (3, 18, 120482232),
    (3, 19, -19609272),
    (3, 21, -23529434),
    (3, 24, -98379162),
    (3, 25, 42147361),
    (4, 9, 272941041),
    (4, 11, 23305245),
    (4, 12, -70035550),
    (4, 14, 241110784),
    (4, 16, 61213470),
    (4, 19, 135393685),
    (4, 21, 30684152),
    (4, 24, 262861878),
    (5, 8, 35191781),
    (5, 10, 321011184),
    (5, 14, -104169070),
    (5, 15, 375111138),
    (5, 20, 418775902),
    (5, 23, -98757666),
    (5, 28, -47386121),
    (6, 7, 381208558),
    (6, 11, -42891154),
    (6, 12, 337627573),
    (6, 13, -146194020),
    (6, 14, -78782821),
    (6, 16, 19526989),
    (6, 19, -86283590),
    (6, 21, 90671780),
    (6, 22, 187815754),
    (6, 23, 202130330),
    (6, 28, -36203985),
    (7, 8, -59663804),
    (7, 9, -100785398),
    (7, 12, 291996527),
    (7, 15, -44507231),
    (7, 16, 55398098),
    (7, 21, 117706337),
    (7, 22, 165147825),
    (8, 11, -56920029),
    (8, 13, 259229426),
    (8, 14, -47902526),
    (8, 15, 57463096),
    (8, 18, 183702360),
    (8, 19, -222235396),
    (8, 20, 41480973),
    (8, 21, 10556935),
    (8, 28, 92833206),
    (9, 14, 379995036),
    (9, 18, -98046773),
    (9, 19, 316558361),
    (9, 28, -53434665),
    (10, 12, -37639343),
    (10, 14, -98004443),
    (10, 15, 294664038),
    (10, 16, -99343891),
    (10, 18, -90692258),
    (10, 19, 60606640),
    (11, 16, 300836730),
    (11, 18, 35163766),
    (11, 28, 237971338),
    (13, 14, -57898358),
    (13, 28, 23099334),
)

def require(condition,message):
    if not condition:
        raise AssertionError(message)

def compact_h_certificate(h_graph,labels=None,entries=H_ORBIT_ENTRIES):
    """Reconstruct the integer H matrix and exact congruence blocks."""
    if labels is None:labels=list(H_LABELS)
    require(labels==list(H_LABELS),"Wrong native H labels")
    n=len(labels)
    require(len(h_graph)==n and all(len(row)==n for row in h_graph),"Wrong graph order")
    require(all(h_graph[i][i]==0 for i in range(n)),"Graph has a loop")
    require(all(h_graph[i][j]==h_graph[j][i] in (0,1) for i in range(n) for j in range(n)),"Graph is not simple symmetric")
    require(sum(map(sum,h_graph))==352,"Wrong H edge count")
    index={v:i for i,v in enumerate(labels)}
    edge=lambda x,y:bool(h_graph[index[x]][index[y]])
    sigma={x:x if x in (0,12,28) else 29-x for x in labels}
    require(all(edge(x,y)==edge(sigma[x],sigma[y]) for x,y in combinations(labels,2)),"Claimed involution is not a graph automorphism")
    M=[[None]*n for _ in labels]
    for x in labels:
        M[index[x]][index[x]]=DIAGONAL
        for y in labels:
            if edge(x,y):M[index[x]][index[y]]=-SCALE
    expected={(x,y) for x,y in combinations(labels,2) if not edge(x,y)
              and (x,y)<=tuple(sorted((sigma[x],sigma[y])))}
    require(len(expected)==91,"Unexpected nonedge orbit count")
    seen=set()
    for x,y,value in entries:
        require((x,y) in expected and (x,y) not in seen,"Wrong or duplicate orbit representative")
        require(type(value) is int,"Nonintegral certificate entry")
        seen.add((x,y))
        for u,v in [(x,y),(y,x),(sigma[x],sigma[y]),(sigma[y],sigma[x])]:
            old=M[index[u]][index[v]]
            require(old is None or old==value,"Inconsistent orbit entries")
            M[index[u]][index[v]]=value
    require(seen==expected and all(x is not None for row in M for x in row),"Incomplete certificate")
    N=lambda x,y:M[index[x]][index[y]]
    require(all(N(x,y)==N(y,x)==N(sigma[x],sigma[y]) for x in labels for y in labels),"Symmetry failure")
    fixed=[0,12,28];reps=[x for x in labels if x<sigma[x]]
    plus=[[2*N(f,g) for g in fixed]+[2*N(f,x) for x in reps] for f in fixed]
    plus += [[2*N(x,f) for f in fixed]+[N(x,u)+N(x,sigma[u]) for u in reps] for x in reps]
    minus=[[N(x,u)-N(x,sigma[u]) for u in reps] for x in reps]
    return {"labels":labels,"matrix":M,"automorphism":sigma,
            "fixed":fixed,"representatives":reps,
            "S":plus,"T":minus,"scale":SCALE,"diagonal":DIAGONAL}

def positive_bareiss_minors(matrix):
    """Leading principal minors by exact fraction-free elimination."""
    M=[row[:] for row in matrix];previous=1;result=[]
    for k in range(len(M)):
        pivot=M[k][k]
        require(pivot>0,"Nonpositive leading principal minor")
        result.append(pivot)
        for i in range(k+1,len(M)):
            for j in range(k+1,len(M)):
                value=M[i][j]*pivot-M[i][k]*M[k][j]
                require(value%previous==0,"Nonintegral Bareiss step")
                M[i][j]=value//previous
        previous=pivot
    return result

def positive_schur_pivots(matrix):
    """Independent rational Schur check on the full 27-by-27 matrix."""
    M=[[Fraction(x) for x in row] for row in matrix];pivots=[]
    for k in range(len(M)):
        pivot=M[k][k];require(pivot>0,"Nonpositive Schur pivot");pivots.append(pivot)
        for i in range(k+1,len(M)):
            for j in range(k+1,len(M)):
                M[i][j]-=M[i][k]*M[k][j]/pivot
    return pivots

def verify_compact_h(h_graph,labels=None,entries=H_ORBIT_ENTRIES):
    result=compact_h_certificate(h_graph,labels,entries)
    result["S_minors"]=positive_bareiss_minors(result["S"])
    result["T_minors"]=positive_bareiss_minors(result["T"])
    result["matrix_Schur_pivots"]=positive_schur_pivots(result["matrix"])
    require(Fraction(DIAGONAL,SCALE)+1==Q,"Wrong q")
    require(Q<Fraction(259828073,50000000),"New bound does not imply previous bound")
    require(Fraction(5)**13*Q**987<L**1000,"First strict power inequality failed")
    require(Fraction(11,2)**13*B**987<L**1000,"Second strict power inequality failed")
    return result

def read_h_from_tex(path):
    source=Path(path).read_text()
    section=source.split(r"\subsubsection{The 27-vertex graph}",1)[1]
    section=section.split(r"\end{longtable}",1)[0]
    labels=list(H_LABELS);index={v:i for i,v in enumerate(labels)}
    adjacency=[[0]*len(labels) for _ in labels];seen=set()
    for left,raw in re.findall(r"^\s*(\d+)&\\\((.*?)\\\)\\\\",section,re.M):
        x=int(left);require(x in labels and x not in seen,"Unexpected graph row");seen.add(x)
        for y in map(int,re.findall(r"\d+",raw)):
            require(x<y and y in labels,"Unexpected graph neighbour")
            adjacency[index[x]][index[y]]=adjacency[index[y]][index[x]]=1
    require(seen==set(labels) and sum(map(sum,adjacency))==352,"Incomplete H graph table")
    return adjacency,labels

if __name__=="__main__":
    adjacency,labels=read_h_from_tex(Path(__file__).resolve().parent.parent/"subsix_supplement.tex")
    result=verify_compact_h(adjacency,labels)
    print(json.dumps({"status":"exact checks passed","orbit_entries":len(H_ORBIT_ENTRIES),
                      "order":len(result["labels"]),"positive_block_orders":[15,12],
                      "positive_matrix_Schur_pivots":len(result["matrix_Schur_pivots"]),
                      "H_matrix":"positive definite",
                      "q":str(Q),"weight":str(WEIGHT),"threshold":str(L),
                      "strict_power_comparisons":[True,True]},indent=2))
