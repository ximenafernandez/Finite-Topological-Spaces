import os
os.system('sage --preparse General.sage')
os.system('mv General.sage.py General.py')
os.system('sage --preparse Homotopy.sage')
os.system('mv Homotopy.sage.py Homotopy.py')

from General import *
from Homotopy import *

#Random spanning tree
def find(C, u):
	if C[u] != u:
		C[u] = find(C, C[u])
	return C[u]

def union(C, u, v):
	u,v = find(C, u), find(C, v)
	C[u] = v

def kruskal(X):
	E = X.cover_relations()
	shuffle(E)
	T = []
	C = {u:u for u in X}
	for e in E:
		if find(C, e[0]) != find(C, e[1]):
			T.append(e)
			union(C, e[0], e[1])
	return T

#Random spanning collapsible
def expand(X,A):
	not_edges = [e for e in X.cover_relations() if not e in A]
	for e in not_edges:
		if is_collapsible(Poset((X.list(), A + [e]))):
			A.append(e)
	return A
			
def spanning_collapsible(X):
	A = kruskal(X)
	while len(expand(X,A)) != len(A):
		A = expand(X,A)
	return A

#Presentation associated to the coloring A
def coloring_presentation(X,A): #A is the list of edges of a spanning collapsible
	gens = [r for r in X.cover_relations() if not r in A]
	F = FreeGroup(len(gens), 'a')
	d = {}
	for i in range(len(gens)):
		d[tuple(gens[i])] = F.gens()[i]
	for j in range(len(A)):
		d[tuple(A[j])] = F.one()

	rels = []
	for x in X.minimal_elements():
		for w in X.maximal_elements():
			l = [u for u in X.upper_covers(x) if u in X.lower_covers(w)]
			ind = -1
			for j in range(len(l)):
				if(not [x,l[j]] in A or not [l[j],w] in A):
					ind = j
		
			if(ind != -1):
				y = l[j]
				for z in l:
					if(z != y):
						s1 = [x,y]
						s2 = [y,w]
						s3 = [x,z]
						s4 = [z,w]
					
						rels.append(d[tuple(s1)] * d[tuple(s2)] * d[tuple(s4)]^-1 * d[tuple(s3)]^-1)
	return  (F / rels)
