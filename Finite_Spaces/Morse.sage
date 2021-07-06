from Finite_Spaces.General       import *
from Finite_Spaces.Homotopy      import *
from Finite_Spaces.Presentations import *

#Morse presentation

def attaching(gens, rels):
	'''
	Computes the original attaching map of each 2-cell in the barycentric subdivision of the standard complex associated to the presentation <gens|rels>.
	INPUT: list of generators and relators in tuple format, e.g x**2 is written as [('x', 2)]
	OUTPUT: dictionary of 2-cells: an ordered list of oriented 1-cells in format (cell, +1 or - 1). 
 	'''
	att = {} #dictionary of attaching maps
	
	#expanded relations of the presentation
	Rels=[]
	for i in range(len(rels)):
		Rels.append([])
		for l in rels[i]:
			if l[1] < 0:
				Rels[i] += [(l[0],-1) for k in range(abs(l[1]))]
			if l[1] > 0:
				Rels[i] += [(l[0],1) for k in range(abs(l[1]))]
 
	for i in range(len(Rels)):
		letters = {}
		R = Rels[i]
		for (l,e) in R:
			letters[l] = 0
		for j in range(len(R)):
			(l,e) = R[j]
			letters[l] += 1
			if(e == 1):
				att[aux_label(i+1,l,'2',letters[l])] = [(l+'2',1),(aux_label(i+1,l,'1',letters[l]),-1),(aux_label(i+1,0,'',j+1),1)]
				att[aux_label(i+1,l,'3',letters[l])] = [(l+'3',1),(aux_label(i+1,0,'',(j+1)%(len(R))+1),-1), (aux_label(i+1,l,'1',letters[l]),1)]
			else:
				att[aux_label(i+1,l,'2',letters[l])] = [(l+'2',-1),(aux_label(i+1,0,'',(j+1)%(len(R))+1),-1),(aux_label(i+1,l,'1',letters[l]),1)]
				att[aux_label(i+1,l,'3',letters[l])] = [(l+'3',-1),(aux_label(i+1,l,'1',letters[l]),-1),(aux_label(i+1,0,'',j+1),1)]
	
	return att

# auxiliary function for attaching_Morse
def write(edge, isolate):
	global new_attaching
	new_attaching = []
	return (aux_write(edge, isolate))

# auxiliary function for attaching_Morse
def aux_write(edge, isolate):
	global new_attaching
	if edge in isolate.keys():
		if isolate[edge] == []:
			return new_attaching
		for edge in isolate[edge]:
			aux_write(edge, isolate)
		return new_attaching
	else:
		new_attaching += [edge]
		return new_attaching

def attaching_Morse(attaching, M, critical_dim_2):
	'''
	Computes the attaching map of each 2-cell in the Morse complex associated to the acyclic matching M.
	INPUT: 
	- attaching: the dictionary of original attaching maps of 2-cells
	- M: acyclic matching 
	- critical_dim_2: the list of critical cells of dimension 2
 	'''
	dim2 = [] # cells of dimension 2 that collapse
	isolate = {}
	for (c1,c2) in M:
		if c2 in attaching.keys(): # c2 of dimension 2
			dim2.append(c2)
			att = [] # new attaching map of c1
			orient = 1
			
			for (cell,e) in attaching[c2]:
				if cell != c1:
					att.append((cell,-e))
				else:
					orient = e
					break
			att = att[::-1] # reverse list
			
			for (cell,e) in reversed(attaching[c2]):
				if cell != c1:
					att.append((cell,-e))
				else:
					break

			if orient == -1:
				att = att[::-1] # reverse list
				for i in range(len(att)): # inverse of every cell
					(cell,e) = att[i]
					att[i] = (cell,-e)
			isolate[(c1,1)] = att
			att_inv = []
			for i in range(len(att)): # inverse of every cell
				(cell,e) = att[i]
				att_inv.append((cell,-e))
				
			att_inv = att_inv[::-1]
			
			isolate[(c1,-1)] = att_inv
		else:
			isolate[(c2,1)] = []
			isolate[(c2,-1)] = []
  
	d = {}
	for c in critical_dim_2:
		rel = []
		for edge in attaching[c]:
			rel += write(edge, isolate)
		d[c] =  rel
	
	return d


def att_to_group(att, gens):
	F = FreeGroup(len(gens), 'a')
	d = {}
	for i in range(len(gens)):
	   d[gens[i]] = F.gens()[i]
	Rels = []
	for c in att.keys():
		aux = F.one()
		for (cell,e) in att[c]:
			aux *= d[cell]^e
		Rels.append(aux)
	return F / Rels


def critical_by_level(X, M): 
	'''
	Return the list of critical points in X associated to the matching M for each level in the face poset
	INPUT:
	- X: poset (face poset of a regular CW)
	- M: (an acyclic) matching 
	'''
	matched = [e[0] for e in M] + [e[1] for e in M]
	edges = X.cover_relations()
	G = DiGraph(edges)
	l = G.level_sets()
	C = []
	for i in range(len(l)):
		cr = []
		for x in l[i]:
			if x not in matched:
				cr.append(x)
		C.append(cr)
	return C

def Morse_presentation(gens, rels, M):
	'''
	Returns the Morse presentation associated to the original presentation <gens|rels> and the matching M.
	'''
	original_attaching = attaching(gens,rels)
	X = presentation_poset(gens,rels)
	critical_cells = critical_by_level(X, M)
	Morse_pres = att_to_group(attaching_Morse(original_attaching, M, critical_cells[2]), critical_cells[1])
	return Morse_pres



#Matchings

def induced_spanning_tree(M, X):
    '''
    Return the induced spanning tree
    INPUT:
    - X face poset of a regular CW-complex
    - M matching in X with unique critical 0-cell
    '''
    T = []
    for e in M:
        P = X.subposet(X.order_ideal([e[1]]))
        if P.height()<2:
            T.append(e)
    return T
	

import random
def greedy_acyclic_matching(X):
	'''
	Greedy algorithm that ouputs a random maximal matching.
	INPUT: X the face poset of regular CW.
	'''
	edges = X.cover_relations()
	in_match = {}
	for v in X.list():
		in_match[v] = False
	seed()
	shuffle(edges)
	M = []
	for e in edges:
		if(in_match[e[0]] or in_match[e[1]]):
			continue
		D = DiGraph(edges)
		D.reverse_edge(e)

		if D.is_directed_acyclic():
			edges.remove(e)
			edges.append([e[1], e[0]])
			M.append(e)
			in_match[e[0]] = True
			in_match[e[1]] = True
	return M

def spanning_matching(X):
	M = []
	n = 0
	while n != 1:
		M = greedy_acyclic_matching(X)
		n = len(critical_by_level(X, M)[0]) == 1
	return M
    
def is_acyclic(X, M):
    edges = X.cover_relations()
    D = DiGraph(edges)
    for e in M:
        D.reverse_edge(e)
    return D.is_directed_acyclic()


#Incidence of critical cells

def critical_points(X, M): #X the face poset of a regular CW, M an acyclic matching
	L = X.list()
	for e in M:
		L.remove(e[0])
		L.remove(e[1])
	return L

#Ouputs the reversed Hasse diagram of X according to M
def reverse(X, M): #X the face poset of a regular CW, M an acyclic matching
	edges = X.cover_relations()
	for e in M:
		edges.remove(list(e))
		edges.append([e[1], e[0]])
	return DiGraph(edges)

#Ouputs the incidence of y in x in K_M
def critical_incidence(gens, rels, M, x, y): #X the face poset of a regular CW (say K), M an acyclic matching, y of height i+1, x of height i. 
	X = presentation_poset(gens,rels)
	incidence = attaching(gens,rels)
	D = reverse(X,M)
	L = D.all_simple_paths(starting_vertices = [x], ending_vertices = [y])
	inc = 0
	for p in L: 
		s = incidence(p[0], p[1])
		for i in range(2, len(p), 2):
			s = s * incidence(p[i], p[i-1]) * incidence(p[i], p[i+1])
			inc = inc + s * (-1)^(len(p)/2 - 1)
	return inc

def critical_CW_incidence(gens, rels, M): #X the face poset of a regular CW, M an acyclic matching
	X = presentation_poset(gens, rels)
	C = critical_by_level(X, M)
	inc = {}
	for i in range(len(C)-1):
		if C[i] != [] and C[i+1] != []:
			for y in C[i+1]:
				for x in C[i]:
					inc[(x, y)] = critical_incidence(gens, rels, M, x, y)
	return inc
