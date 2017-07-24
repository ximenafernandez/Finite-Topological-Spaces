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
            union(C, e[0],e[1])
    return T

#Random spanning collapsible
def expand(X,A):
    not_edges=[e for e in X.cover_relations() if not e in A]
    for e in not_edges:
        if is_collapsible(Poset((X.list(), A+[e]))):
            A.append(e)
    return A
            
def spanning_collapsible(X):
    A=kruskal(X)
    while len(expand(X,A))!=len(A):
        A=expand(X,A)
    return A

#Presentation associated to the coloring A
def coloring_presentation(X,A): #A is the list of edges of a spanning collapsible
    gens=[r for r in X.cover_relations() if not r in T]
    F=FreeGroup(len(gens), 'a')
    dict={}
    for i in range(len(gens)):
        dict[tuple(gens[i])]=F.gens()[i]
    for j in range(len(T)):
        dict[tuple(T[j])]=F.one()
    
    rels=[]
    for x in X.minimal_elements():
        l=X.upper_covers(x)
        for y,z in Subsets(l,2):
            for w in X.maximal_elements():
                if X.is_gequal(w,y) and X.is_gequal(w,z):
                    s1=[x,y]
                    s2=[y,w]
                    s3=[x,z]
                    s4=[z,w]
                    rels.append(dict[tuple(s1)]*dict[tuple(s2)]*dict[tuple(s4)]^-1*dict[tuple(s3)]^-1)             
    return  F/rels
