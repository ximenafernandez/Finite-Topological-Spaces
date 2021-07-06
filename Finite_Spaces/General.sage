#Opposite of a finite space
def op(X):
	elms = X.list()
	rels = [r[::-1 ] for r in X.cover_relations()]
	return Poset((elms,rels))

#Join
def join (X, Y): 
	return X.ordinal_sum(Y)

#Cartesian Product
def cartesian_product (X, Y):
	return X.product(Y)

#Quotient
def quotient_poset(X, A): #X poset, A the list of elements of the subposet
	open_A = X.order_ideal(A)
	closed_A = X.order_filter(A)
	if set(open_A).intersection(closed_A) != set(A):
		return 'ERROR: X/A is not T_0'
	notinA = [x for x in X.list() if not x in A]
	elms = notinA + [A[0]]
	rels=[[x, y] for x in notinA for y in notinA if X.covers(x, y)] + [[A[0], x] for a in A for x in notinA if X.covers(a, x)] + [[x, A[0]] for a in A for x in notinA if X.covers(x, a)]
	return Poset((elms, rels))

#Wedge
def wedge(X, Y, x, y):
	return quotient_poset(X.disjoint_union(Y), [(0, x), (1, y)])

