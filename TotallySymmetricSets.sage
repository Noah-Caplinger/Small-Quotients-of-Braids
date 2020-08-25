
### I'm re-doing some of the stuff in TSSinPSL27.sage
### because it was starting to get messy.

import sage.groups.perm_gps.permgroup

########################################################### 
##                    Basic Groups
###########################################################

PSL27 = PermutationGroup([[(1,2,3),(5,6,7)],[(2,4),(3,5)]])
PSL28 = PermutationGroup([[(1,2),(3,4),(6,7),(8,9)],[(1,3,2),(4,5,6),(7,8,9)]])
PSL211 = PermutationGroup([[(2,10),(3,4),(5,9),(6,7)],[(1,2,11),(3,5,10),(6,8,9)]])
PSL213 = PermutationGroup([[(1,12),(2,6),(3,4),(7,11),(9,10),(13,14)],[(1,6,11),(2,4,5),(7,8,10),(12,14,13)]])
PSL217 = PermutationGroup([[(1,16),(2,8),(3,11),(5,10),(6,14),(7,12),(9,15),(17,18)],[(1,8,15),(2,11,7),(3,4,10),(5,14,9),(6,12,13),(16,18,17)]])
PSL219 = PermutationGroup([[(1,18),(2,9),(3,6),(4,14),(5,15),(7,8),(10,17),(11,12),(13,16),(19,20)],[(1,9,17),(2,6,8),(3,14,5),(4,15,13),(10,12,16),(18,20,19)]])
PSL216 = PermutationGroup([[(1,2),(3,4),(5,6),(7,9),(10,11),(12,13),(14,16),(15,17)],[(1,3,2),(4,5,7),(6,8,10),(11,12,14),(13,15,17)]])
PSL33 = PermutationGroup([[(2,4),(3,5),(6,8),(10,11)],[(1,2,3),(5,6,7),(8,9,10),(11,12,13)]])
PSL34 = PermutationGroup([[(1,2),(4,6),(5,7),(8,12),(9,14),(10,15),(11,17),(13,19)],[(2,3,5,4),(6,8,13,9),(7,10,16,11),(12,18),(14,20,21,15),(17,19)]])
PSL227 = PermutationGroup([[(1,2),(3,4),(5,7),(6,8),(9,12),(10,13),(11,15),(14,19),(16,20),(17,21),(18,22),(23,25),(24,26),(27,28)],[(2,3,5),(4,6,9),(7,10,14),(8,11,16),(12,17,20),(13,18,23),(15,19,22),(21,24,27),(25,26,28)]])

G22_prime = PermutationGroup([[(2,3),(4,6),(5,8),(7,11),(9,13),(10,15),(12,14),(16,20),(17,22),(18,23),(24,27),(25,28)],[(1,2,4,7,12,17),(3,5,9,14,19,22),(6,10,13,18,24,23),(8,11,16,21,26,28),(20,25,27)]])

PGL29 = PermutationGroup([[(1,2),(3,5),(4,6),(7,8),(9,10)],[(2,3,4),(5,7,8),(6,9,10)]])

S2 = SymmetricGroup(2)
S3 = SymmetricGroup(3)
S4 = SymmetricGroup(4)
S5 = SymmetricGroup(5)
S6 = SymmetricGroup(6)
S7 = SymmetricGroup(7)
S8 = SymmetricGroup(8)

A3 = AlternatingGroup(3)
A4 = AlternatingGroup(4)
A5 = AlternatingGroup(5)
A6 = AlternatingGroup(6)
A7 = AlternatingGroup(7)
A8 = AlternatingGroup(8)

#http://sci.tech-archive.net/Archive/sci.math/2006-12/msg07456.html
M_10 = PermutationGroup([[(2,3,4),(5,6,7),(8,9,10)],[(2,5,8),(3,6,9),(4,7,10)],[(3,5,4,8),(6,7,10,9)],[(3,6,4,10),(5,9,8,7)],[(1,2),(5,8),(6,7),(9,10)]]) 
M_11 = PermutationGroup([[(2,10),(4,11),(5,7),(8,9)],[(1,4,3,8),(2,5,6,9)]])

########################################################### 
##                    Basic operations
###########################################################

def commutes_with_list(g, alist):
	for h in alist:
		if g*h != h*g:
			return False
	return True

def braids_with(a,b):
	return a*b*a == b*a*b

def commutes_with(a,b):
	return a*b == b*a


# takes in a group and alist, and returns a list containing one 
# representative from each conjugacy class of lists appearing in alist.

def up_to_conjugation(G,alist):
	toReturn = []
	while alist != []:
		S = alist[0]
		toReturn.append(S)
		for g in G:  #G, not the conjugacy class of x[0], because I think that requires the same overhead
			x = [g*S[i]/g for i in range(0,len(S))]
			if x in alist:
				alist.remove(x)
	return toReturn


########################################################### 
##         Check if all permutations are realized
###########################################################

## transpositions of the form (1,x)
transpositions_1_x = []

for i in range(2, 25):
	transpositions_1_x.append(PermutationGroupElement([(1,i)]))

## returns a list of generating transpositions of S (1 2), (1 3), (1 4), etc 
## So [1,2,3] returns [[2,1,3], [3,2,1]], ie the permutations (1 2) and (1 3)

def get_transpositions(S):
	toReturn = []
	relavent_generating_transpositions = transpositions_1_x[:len(S) - 1]
	for transposition in relavent_generating_transpositions:
		x = [S[transposition(i + 1) - 1] for i in range(len(S))]
		toReturn.append(x)
	return toReturn


# G is the ambient group, S is the set in question
# NOTE: this only checks that every permutation of the set is realized via conjugation
# and NOT the commutativity of the elements of S

def all_permutations_realized(G, S):

	generating_transpositions = get_transpositions(S)

	for g in G:
		x = [g*S[i]/g for i in range(len(S))]
		if x in generating_transpositions:
			generating_transpositions.remove(x)

	return generating_transpositions == []



########################################################### 
##                   Finding TSS
###########################################################


## Finds all totally symmetric sets in G of size 2

def TSS_size_2(G):
	toReturn = []
	conj_classes = G.conjugacy_classes()
	for c in conj_classes:
		c = c.list()
		for i in range(0,len(c)):
			for j in range(i+1,len(c)):
				if c[i]*c[j] != c[j]*c[i]:
					continue
				TS_canidate = [c[i],c[j]]
				if all_permutations_realized(G,TS_canidate):
					toReturn.append(TS_canidate)
	return toReturn

#sometimes faster; does not respect ordering; not used

def alternate_TSS_size_2(G):
	withRepeats = []
	for first_elt in G:
		for conj_elt in G:
			second_elt = conj_elt*first_elt*conj_elt.inverse()
			if second_elt == first_elt:
				continue
			if not commutes_with(first_elt, second_elt):
				continue
			if conj_elt*second_elt*conj_elt.inverse() == first_elt: ##TSS condition
				withRepeats.append([first_elt,second_elt])
	##delete repeats
	noRepeats = []
	for TSS in withRepeats:
		if TSS not in noRepeats and [TSS[1],TSS[0]] not in noRepeats:
			noRepeats.append(TSS)
	return noRepeats






## Constructs TSS recursively

def TSS_size_n(G,n):
	if n == 2:
		#return up_to_conjugation(G,all_TSS_size_2(G))
		return TSS_size_2(G)

		## this is bad, because it doesnt respect ordering of G
		## which is assumed in "index shenanigans"
		#return alternate_TSS_size_2(G) 


	smaller_TSSs = TSS_size_n(G,n-1)

	toReturn = []

	## now to construct a larger TSS

	for small_TSS in smaller_TSSs:
		last_elt = small_TSS[-1]
		conj_class_of_last_elt = G.conjugacy_class(last_elt).list()

		## indexing shenanigans to inherit ordering from the conjugacy class
		for i in range(conj_class_of_last_elt.index(last_elt) + 1, len(conj_class_of_last_elt)):
			g = conj_class_of_last_elt[i]
			if not commutes_with_list(g,small_TSS) or g in small_TSS:
				continue
			if all_permutations_realized(G,small_TSS + [g]):
				toReturn.append(small_TSS + [g])

	#return up_to_conjugation(toReturn)
	return toReturn


## We might also know the result of the inductive step already:

def TSS_size_n_smaller_known(G,n,smaller_TSSs):
	if n == 2:
		#return up_to_conjugation(G,all_TSS_size_2(G))
		return TSS_size_2(G)

	toReturn = []

	## now to construct a larger TSS

	for small_TSS in smaller_TSSs:
		last_elt = small_TSS[-1]
		conj_class_of_last_elt = G.conjugacy_class(last_elt).list()

		## indexing shenanigans to inherit ordering from the conjugacy class
		for i in range(conj_class_of_last_elt.index(last_elt) + 1, len(conj_class_of_last_elt)):
			g = conj_class_of_last_elt[i]
			if not commutes_with_list(g,small_TSS) or g in small_TSS:
				continue
			if all_permutations_realized(G,small_TSS + [g]):
				toReturn.append(small_TSS + [g])

	#return up_to_conjugation(toReturn)
	return toReturn


########################################################### 
##               Putting things together
###########################################################

## returns a list of list of TSS, where the 
## 0th position has the TSS of size 2, 1st has size 3, etc

def TSS_all_sizes(G):
	toReturn = [TSS_size_2(G)]
	size = 2
	while toReturn[-1] != []:
		size += 1
		toReturn.append(TSS_size_n_smaller_known(G,size,toReturn[-1]))
	for i in range(len(toReturn)):
	 	toReturn[i] = up_to_conjugation(G,toReturn[i])
	return toReturn


## find order set, given all sizes of TSS:
def order_set(all_sizes):
	toReturn = set()
	for size in range(len(all_sizes)-1):
		for TSS in all_sizes[size]:
			##find the order of TSS[0]
			P = PermutationGroup([TSS[0].cycle_tuples()])
			toReturn.add(P.order())
	return toReturn



def print_TSS_info(G):
	all_sizes = TSS_all_sizes(G)

	for i in range(0,len(all_sizes) - 1):
		print("All TSS of size " + str(i+2) + ": ")
		print(all_sizes[i])
		print("")	

	print("Max size: " + str(len(all_sizes)))
	print("")

	print("Order set: ", order_set(all_sizes))










