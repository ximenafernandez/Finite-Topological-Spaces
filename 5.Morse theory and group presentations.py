
new_adj = []

def write(edge, isolate):
    global new_adj
    #print edge, new_adj
    if edge in isolate.keys():
        if isolate[edge] == []:
            return new_adj
        for edge in isolate[arista]:
            write(edge, isolate)
        return new_adj
    else:
        new_adj += [edge]
        return new_adj

def adj_Morse(adj, matching, critics):
    global new_adj
    dim2 = [] # cells of dimension 2 that collapses
    isolate={}
    for (c1,c2) in matching:
        if c2 in adj.keys(): # c2 dimension 2
            dim2.append(c2)
            att = [] # new writting of c1
            orient = 1
            
            #print (c1,c2)
            
            for (cel,e) in adj[c2]:
                if cel != c1:
                    att.append((cel,-e))
                else:
                    orient = e
                    break
            att = att[::-1] # revert list
            
            #print att
            
            for (cel,e) in reversed(adj[c2]):
                if cel != c1:
                    att.append((cel,-e))
                else:
                    break
                    
            #print att
            
            if orient == -1:
                att = att[::-1] # revert list
                for i in range(len(att)): # invert every cell
                    (cel,e) = att[i]
                    att[i] = (cel,-e)
            isolate[(c1,1)] = att
            att_inv = []
            for i in range(len(att)): # invert every cell
                (cel,e) = att[i]
                att_inv.append((cel,-e))
                
            att_inv = att_inv[::-1]
            
            isolate[(c1,-1)] = att_inv
        else:
            isolate[(c2,1)] = []
            isolate[(c2,-1)] = []

    dict = {}
    for c in critics:
        rel = []
        for edge in adj[c]:
            new_adj = []
            rel += write(edge, isolate)
        dict[c] =  rel
    
    return dict

#Greedy algorithm that ouputs a random maximal matching

import random
def greedy_acyclic_matching(X): #X the face poset of regular CW.
    edges=X.cover_relations()
    available_edges=X.cover_relations()
    seed()
    shuffle(available_edges)
    M=[]
    while available_edges!=[]:
        n=0
        for e in available_edges:
            D=DiGraph(edges)
            D.reverse_edge(e)
            if D.is_directed_acyclic():
                n=1
                edges.remove(e)
                edges.append([e[1],e[0]])
                M.append(e)
                available_edges.remove(e)
                l=[]
                for r in available_edges:
                    if r[0]==e[0] or r[1]==e[0] or r[0]==e[1] or r[1]==e[1]:
                        l.append(r)
                for r in l:
                    available_edges.remove(r)
        if n==0:return M            
    return M

def critical_points(X,M): #X the face poset of a regular CW, M an acyclic matching
    L=X.list()
    for e in M:
        L.remove(e[0])
        L.remove(e[1])
    return L

#Ouputs the reversed Hasse diagram of X according to M
def reverse(X,M): #X the face poset of a regular CW, M an acyclic matching
    edges=X.cover_relations()
    for e in M:
        print e
        edges.remove(list(e))
        edges.append([e[1],e[0]])
    return DiGraph(edges)

def criticals_by_level(X,M): #X the face poset of a regular CW, M an acyclic matching
    matched=[e[0] for e in M]+[e[1] for e in M]
    edges=X.cover_relations()
    G=DiGraph(edges)
    l=G.level_sets()
    C=[]
    for i in range(len(l)):
        c=[]
        for x in l[i]:
            if x not in matched:
                c.append(x)
        C.append(c)
    return C

#Ouputs the incidence of y in x in K_M
def critical_incidence(gens,rels,M,x,y): #X the face poset of a regular CW (say K), M an acyclic matching, y of height i+1, x of height i. 
	X=pres_poset(gens,rels)
	incidence=incidence(gens,rels)
    D=reverse(X,M)
    L=D.all_simple_paths(starting_vertices=[x], ending_vertices=[y])
    inc=0
    for p in L: 
		s=incidence(p[0],p[1])
		for i in range(2,len(p),2):
		s=s*incidence(p[i],p[i-1])*incidence(p[i],p[i+1])
        #print len(p), len(p)/2-1
        inc=inc+s*(-1)^(len(p)/2-1)
    return inc


def critical_CW_incidence(gens,rels,M): #X the face poset of a regular CW, M an acyclic matching
	X=pres_poset(gens,rels)
    C=criticals_by_level(X,M)
    inc={}
    for i in range(len(C)-1):
        if C[i]!=[] and C[i+1]!=[]:
            for y in C[i+1]:
                for x in C[i]:
                    inc[(x,y)]= incidence(gens,rels,M,x,y)
    return inc

def adj_to_group(adj,gens):
    F=FreeGroup(len(gens), 'a')
    dict={}
    for i in range(len(gens)):
       dict[gens[i]]=F.gens()[i]
    Rels=[]
    for c in adj.keys():
        aux = F.one()
        for (celda,e) in adj[c]:
            aux *= dict[celda]^e
        Rels.append(aux)
    return F/Rels

def presentation_matching(X):
    M=[]
    n=0
    while n!=1:
        M=greedy_acyclic_matching(X)
        n=len(criticals_by_level(X,M)[0])==1
    return M

def len_rels(G):
    return sum(sum(abs(e) for (l,e) in r.syllables())for r in G.relations())
