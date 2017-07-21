#Osaki reductions

def is_down_O_reduction(X,x):
    if len(X.order_ideal([x]))==1: return False
    for y in X:
        intersection=Set(X.order_ideal([x])).intersection(Set(X.order_ideal([y]))) 
        if not(is_contractible(X.subposet(list(intersection))) or len(intersection)==0): return False
    return True

def is_up_O_reduction(X,x):
    if len(X.order_filter([x]))==1: return False
    for y in X:
        intersection=Set(X.order_filter([x])).intersection(Set(X.order_filter([y]))) 
        if not(is_contractible(X.subposet(list(intersection))) or len(intersection)==0): return False
    return True

def down_O_core(X):
    if len(X.list())==1: return X
    for x in X.list():
        if is_down_O_reduction(X,x):
            X=quotient(X,X.order_ideal([x]))
            return down_O_core(X)
    return X

def up_O_core(X):
    return down_O_core(op(X))

#Middle reduction 

def is_middle_reduction(X,a,b):
    lower_intersection=Set(X.order_ideal([a])).intersection(Set(X.order_ideal([b])))
    if len(lower_intersection)!=1:
        return False
    Fa_minus_Fb=Set(X.order_filter([a])).difference (Set(X.order_filter([b])))
    for x in Fa_minus_Fb:
		if len(Set(X.order_ideal([b])).intersection(Set(X.order_ideal([x]))))!=1:
            return False
    Fb_minus_Fa=Set(X.order_filter([b])).difference (Set(X.order_filter([a])))
    for x in Fb_minus_Fa:
		if len(Set(X.order_ideal([a])).intersection(Set(X.order_ideal([x]))))!=1:
            return False
    return True
    
def is_middle_op_reduction(X,a,b):
    return is_middle_reduction(op(X),a,b)    

def middle_core(X):
    height1=[x for x in X.list() if not x in X.maximal_elements() and not x in X.minimal_elements()]
    for S in Set(height1).subsets(2):
        l=list(S)
        if is_middle_reduction(X,l[0], l[1]) or is_middle_op_reduction(l[0], l[1], X):
            X=quotient(X,[l[0],l[1]])
            return middle_core(X)
    return X

#Edge reduction

def is_down_reducible_edge(X,e): #e=[e[0], e[1]]
    if e[1] not in X.maximal_elements(): return False
    Y=remove_edge(U(X, e[1]), e)
    return is_contractible(Y)	
    
def is_up_reducible_edge(X,e):
    if e[0] not in X.minimal_elements(): return False
    Y=remove_edge(F(X, e[0]), e)
    return is_contractible(Y)

def is_reducible_edge(X,e):
    return is_up_reducible_edge(X,e) or is_down_reducible_edge(X,e)

def edge_core(X):
    for e in X.cover_relations():
        if is_reducible_edge(X,e):
			X=remove_edge(X,e)
            return edge_core(X)
    return X

#Random Reduction

from random import shuffle

def random_core(X):
    if len(X.list())==1:
        print 'se 3-deforma a un punto'
        return X
    count=range(7)
    shuffle(count)
    for i in range(7):
        Y=random_reduction(X,count[i])
        if X!=Y:
            return random_core(Y)
    return X



def random_reduction(X,j):

    elms=X.list()
    
    if j==0:
        shuffle(elms)
        for x in elms:
            if is_weak_point(X,x):
                print x, 'weak point'
                Y=remove_point(X,x)
                return Y
        return X
        
    if j==1:
        M=X.maximal_elements()
        shuffle(M)
        for S in Set(M).subsets(2):
            l=list(S)
            if is_qc_reduction(X,l[0],l[1]):
                print l[0],l[1], 'qc-reduction'
                Y=quotient(X,[l[1],l[0]])
                return Y
        return X

    if j==2:
        m=X.minimal_elements()
        shuffle(m)
        for S in Set(m).subsets(2):
            l=list(S)
            if is_qc_op_reduction(X,l[0],l[1]):
                print l[0],l[1], 'qc-op-reduction'
                Y=quotient(X,[l[1],l[0]])
                return Y
        return X

    if j==3:
        mid=[x for x in X.list() if not (x in X.maximal_elements() or x in X.minimal_elements())]
        shuffle(mid)
        for S in Set(mid).subsets(2):
            l=list(S)
            if is_middle_reduction(X,l[0],l[1]):
                print l[0],l[1], 'middle-reduction'
                Y=quotient(X,[l[1],l[0]])
                return Y
        return X


    if j==4:
        E=X.cover_relations()
        shuffle(E)
        for e in E:
            if is_reducible_edge(X,e):
                print e, 'edge-reduction'
                Y=remove_edge(X,e)
                return Y
        return X

    if j==5:
        shuffle(elms)
        for x in elms:
            if is_down_O_reduction(X,x):
                print x, 'down_o_reduction'
                Y=quotient(X,X.order_ideal([x]))
                return Y
        return X
     
    if j==6:
        shuffle(elms)
        for x in elms:
            if is_up_O_reduction(X,x):
                print x, 'up_o_reduction'
                Y=quotient(X,X.order_filter([x]))
                return Y
        return X

