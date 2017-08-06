import os
os.system('sage --preparse General.sage')
os.system('mv General.sage.py General.py')
os.system('sage --preparse Homotopy.sage')
os.system('mv Homotopy.sage.py Homotopy.py')
os.system('sage --preparse Presentations.sage')
os.system('mv Presentations.sage.py Presentations.py')

from General import *
from Homotopy import *
from Presentations import *

#Morse 

def attaching(gens, rels):
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

def write(edge, isolate):
	global new_attaching
	new_attaching = []
	return (aux_write(edge, isolate))

def aux_write(edge, isolate):
	global new_attaching
	if edge in isolate.keys():
		if isolate[edge] == []:
			return new_attaching
		for edge in isolate[edge]:
			write(edge, isolate)
		return new_attaching
	else:
		new_attaching += [edge]
		return new_attaching

def attaching_Morse(attaching, matching, critics_dim_2):
	dim2 = [] # cells of dimension 2 that collapse
	isolate = {}
	for (c1,c2) in matching:
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
	for c in critics_dim_2:
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

def critical_by_level(X, M): #X the face poset of a regular CW, M an acyclic matching
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

def morse_presentation(gens, rels, M):
	original_attaching = attaching(gens,rels)
	X = presentation_poset(gens,rels)
	criticals = critical_by_level(X, M)
	morse_attaching = attaching_Morse(original_attaching, M, criticals[2])
	return att_to_group(morse_attaching, criticals[1])

def total_relator_len_(G):
	return sum(sum(abs(e) for (l,e) in r.syllables()) for r in G.relations())

#Matchings

#Greedy algorithm that ouputs a random maximal matching
import random
def greedy_acyclic_matching(X): #X the face poset of regular CW.
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
		print e
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
