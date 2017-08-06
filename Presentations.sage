#Poset associated to a group presentation

def aux_label(a, b, ind, c):
	return ('c' + str(a) + '_' + str(b) + ind + '_' + str(c))

def presentation_poset(gens,rels):
   
	#expanded relations
	Rels = []
	for i in range(len(rels)):
		Rels.append([])
		for l in rels[i]:
			if l[1] < 0:
				Rels[i] += [(l[0], -1) for k in range(abs(l[1]))]
			if l[1] > 0:
				Rels[i] += [(l[0], 1) for k in range(abs(l[1]))]
	
	#V is the set of elements of the presentation poset
	V = range(len(rels)+1) + [x + str(i) for x in gens for i in range(1,4)]
	
	#E is the list of cover relations of the presentaton poset
	E = []
	#edges associated to the generators
	for x in gens:
		E += [[0, x +'2'], [0, x +'3'], [x +'1', x +'2'], [x +'1', x +'3']]
	
	for i in range(len(Rels)):
		letters = {}
		for x in gens: letters[x] = 0
		for j in range(len(Rels[i])):
			letters[Rels[i][j][0]] = letters[Rels[i][j][0]] + abs(Rels[i][j][1])
		total_letters = sum([letters[x] for x in gens])

		V = V + [aux_label(i+1, x, str(j), k) for x in gens if letters[x] != 0 for j in range(1, 4) for k in range(1, letters[x] + 1)] + [aux_label(i+1, 0, '', k) for k in range(1, total_letters + 1)]

		#the model of D^2 i associated to the cell corresponding to the relator i:
		
		#edges between the indicator i of the relator and the minimals of the cycle
		
		E = E + [[i+1, aux_label(i+1, 0, '', j+1)] for j in range(total_letters)] + [[i+1, aux_label(i+1, x, '1', j+1)] for x in gens for j in range(letters[x])]

		# edges between the models of S^1 associated to generators and relators
		
		for x in gens:
			for k in range(1,4):
				for j in range(letters[x]):
					E.append([x + str(k), aux_label(i+1, x, str(k), j+1)])
		for j in range(len(Rels[i])):
			E.append([0, aux_label(i+1, 0, '', j+1)])
			#edges of the cycle which do not start in 0
		for l in Rels[i]:
			for j in range(letters[l[0]]):
				E += [[aux_label(i+1, l[0], '1', j+1), aux_label(i+1, l[0],'2', j+1)],[aux_label(i+1, l[0], '1', j+1),aux_label(i+1, l[0], '3', j+1)]]
			#edges of the cycle starting in 0
		cont = {}
		for x in gens: cont[x] = 0
		
		for j in range(len(Rels[i])):
			cont[Rels[i][j][0]] += 1
			
			#edges to the 'right'
			if Rels[i][j][1] > 0:
				E.append([aux_label(i+1, 0, '', sum([cont[x] for x in gens])), aux_label(i+1, Rels[i][j][0], '2', cont[Rels[i][j][0]])])
			else:
				E.append([aux_label(i+1, 0, '', sum([cont[x] for x in gens])), aux_label(i+1, Rels[i][j][0], '3', cont[Rels[i][j][0]])])
			#edges to the 'left'
			if j != 0:
				if Rels[i][j-1][0] == Rels[i][j][0]:
					if Rels[i][j-1][1] > 0:
						E.append([aux_label(i+1, 0, '', sum([cont[x] for x in gens])), aux_label(i+1, Rels[i][j-1][0], '3', cont[Rels[i][j-1][0]]-1)])
					else:
						E.append([aux_label(i+1, 0, '', sum([cont[x] for x in gens])), aux_label(i+1,Rels[i][j-1][0], '2', cont[Rels[i][j-1][0]]-1) ]) 
				else:
					if Rels[i][j-1][1] > 0:
						E.append([aux_label(i+1, 0, '', sum([cont[x] for x in gens])), aux_label(i+1,Rels[i][j-1][0], '3', cont[Rels[i][j-1][0]])])
					else:
						E.append([aux_label(i+1, 0, '', sum([cont[x] for x in gens])), aux_label(i+1,Rels[i][j-1][0], '2', cont[Rels[i][j-1][0]]) ])
		
		#j=0
		n = len(Rels[i]) - 1
		if Rels[i][n][1] > 0:
			E.append([aux_label(i+1, 0, '', 1), aux_label(i+1, Rels[i][n][0], '3', cont[Rels[i][n][0]])])
		else:
			E.append([aux_label(i+1, 0, '', 1), aux_label(i+1, Rels[i][n][0], '2', cont[Rels[i][n][0]])])

	return Poset((V, E))
