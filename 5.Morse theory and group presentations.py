#Morse

def attaching(gens,rels):
    att={} #dictionary of attaching maps
    
    #relaciones de la presentación expandidas
    Rels=[]
    for i in range(len(rels)):
        Rels.append([])
        for l in rels[i]:
            if l[1]<0:
                Rels[i]=Rels[i]+[(l[0],-1) for k in range(abs(l[1]))]
            if l[1]>0:
                Rels[i]=Rels[i]+[(l[0],1) for k in range(abs(l[1]))]
 
    for i in range(len(Rels)):
        letters={}
        R = Rels[i]
        for (l,e) in R:
            letters[l]=0
        for j in range(len(R)):
            (l,e) = R[j]
            letters[l] += 1
            if(e == 1):
                att['c'+str(i+1)+'_'+l+'2'+'_'+str(letters[l])] = [(l+'2',1),('c'+str(i+1)+'_'+l+'1'+'_'+str(letters[l]),-1),('c'+str(i+1)+'_'+'0'+'_'+str(j+1),1)]
                att['c'+str(i+1)+'_'+l+'3'+'_'+str(letters[l])] = [(l+'3',1),('c'+str(i+1)+'_'+'0'+'_'+str((j+1)%(len(R))+1),-1), ('c'+str(i+1)+'_'+l+'1'+'_'+str(letters[l]),1)]
            else:
                att['c'+str(i+1)+'_'+l+'2'+'_'+str(letters[l])] = [(l+'2',-1),('c'+str(i+1)+'_'+'0'+'_'+str((j+1)%(len(R))+1),-1),('c'+str(i+1)+'_'+l+'1'+'_'+str(letters[l]),1)]
                att['c'+str(i+1)+'_'+l+'3'+'_'+str(letters[l])] = [(l+'3',-1),('c'+str(i+1)+'_'+l+'1'+'_'+str(letters[l]),-1),('c'+str(i+1)+'_'+'0'+'_'+str(j+1),1)]
    
    return att


new_attaching = []

def write(edge, isolate):
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
    global new_attaching
    dim2 = [] # cells of dimension 2 that collapse
    isolate={}
    for (c1,c2) in matching:
        if c2 in attaching.keys(): # c2 of dimension 2
            dim2.append(c2)
            att = [] # new attaching map of c1
            orient = 1
            
            #print (c1,c2)
            
            for (cell,e) in attaching[c2]:
                if cell != c1:
                    att.append((cell,-e))
                else:
                    orient = e
                    break
            att = att[::-1] # revertir lista
            
            #print att
            
            for (cell,e) in reversed(attaching[c2]):
                if cell != c1:
                    att.append((cell,-e))
                else:
                    break
                    
            #print att
            
            if orient == -1:
                att = att[::-1] # reverse
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
  
    dict = {}
    for c in critics_dim_2:
        rel = []
        for edge in attaching[c]:
            new_attaching = []
            rel += write(edge, isolate)
        dict[c] =  rel
    
    return dict

def att_to_group(att,gens):
    F=FreeGroup(len(gens), 'a')
    dict={}
    for i in range(len(gens)):
       dict[gens[i]]=F.gens()[i]
    Rels=[]
    for c in att.keys():
        aux = F.one()
        for (cell,e) in att[c]:
            aux *= dict[cell]^e
        Rels.append(aux)
    return F/Rels

def critical_by_level(X,M): #X the face poset of a regular CW, M an acyclic matching
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

def morse_presentation(gens, rels, M):
	original_attaching=attaching(gens,rels)
	X=presentation_poset(gens,rels)
	critical_dim_2=critical_by_level(X,M)[2]
	critical_dim_1=critical_by_level(X,M)[1]
	new_attaching=attaching_Morse(original_attaching, M, critical_dim_2)
	return att_to_group(new_attaching, critical_dim_1)

def total_relator_len_(G):
    return sum(sum(abs(e) for (l,e) in r.syllables())for r in G.relations())

####################################################################
#Matchings

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

def matching_spanning(X):
    M=[]
    n=0
    while n!=1:
        M=greedy_acyclic_matching(X)
        n=len(critical_by_level(X,M)[0])==1
    return M

############################################################


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






