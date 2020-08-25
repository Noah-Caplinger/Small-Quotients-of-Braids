
# This file is stores functions necessary to check the existence
# of certian homomorphisms from Braid groups

########################################################### 
##                   Sage Import Shenanigans
###########################################################

import os
os.system('sage --preparse TotallySymmetricSets.sage')
os.system('mv TotallySymmetricSets.sage.py TotallySymmetricSets.py')
from TotallySymmetricSets import *


########################################################### 
##                   Maps from B_6'
###########################################################

## This function returns False if we're sure there can't be a 
## homomorphism B_6' -> G, and True otherwise. The basic strategy
## is to try to find consistent images of generators of B_6'
## if we don't, there can't be a homomorphism.

def possible_homomorphism_from_com_B6(G):

	## We will show that the generators of the Lin presentation cannot be sent to 
	## images in PSL(2,7) which are coherent with the relations in B_6'

	## From TSS theory, we know {c1,c3} must be sent to another
	## TSS of size 2 in PSL(2,7).

	TSS_in_G_size_2 = up_to_conjugation(G,TSS_size_2(G))

	possible_homomorphism = False

	for TSS in TSS_in_G_size_2:
		## Because any permutation of TSS can be represented by conjugation
		## by an element, we may assume without loss of generality that 
		## that the c1 is sent to TSS[0], and c2 to TSS[1].
		## If not, simply conjugate the map so that it does.

		c1 = TSS[0]
		c3 = TSS[1]

		canidates_c2= []

		## we will now search for possible images of c2 in PSL(2,7)
		## any canidate must be in the conjugacy class of the image of c1
		## and must braid with the images of c1 and c3

		for c2 in G.conjugacy_class(TSS[0]).list():
			if braids_with(c2,TSS[0]) and braids_with(c2,TSS[1]):
				canidates_c2.append(c2)

		## we've now found all canidates for c2, and will start looking 
		## for potiental images of u = s2s1^-1. These elements must be in the 
		## conjugacy class of the image of s4s3^-1 = (s4s1^-1) (s3s1^-1)^-1, 
		## (ie c2 * c1^-1), 

		for c2 in canidates_c2:
			canidates_u= []
			for u in G.conjugacy_class(c2 * c1.inverse()).list():

				## the following are from the Lin relations
				w = u * c1 * u.inverse()     ## 1st relation
				v  = c2.inverse() * u * c2   ## 5th relation (a)

				r2 = (u*w*u.inverse() == w*w*c1.inverse()*w)
				r3 = (v*c1*v.inverse() == c1.inverse()*w)
				r4 = (v*w*v.inverse() == c1.inverse()*w*c1.inverse()*w*c1.inverse()*w*c1.inverse()*c1.inverse()*w)
				r5b = (u*c3 == c3*v)
				r6a = (v*c2 == c2*u.inverse()*v)
				r6b = (v*c3 == c3*u.inverse()*v)

				if r2 and r3 and r4 and r5b and r6a and r6b:
					canidates_u.append(u)

			if canidates_u != []:
				possible_homomorphism = true
	return possible_homomorphism



def print_homomorphism_from_com_B6(G, name):
	possible_homomorphism = possible_homomorphism_from_com_B6(G)

	if not possible_homomorphism:
		print("\nThere is no non-trivial homomorphism B_6' -> " + name + "\n")
	else:
		print("\nThere might be a non-trivial homomorphism B_6' -> " + name + "\n")

########################################################### 
##                   Maps from B_6
###########################################################

## Although we are considering a homomorphism from B_6 instead of 
## B_6', the basic idea will be the same. We know that the image of 
## {s1, s3, s5} is totally symmetric of size 3 in M_10, and we will 
## try to find possible images for s2 and s4.



def possible_homomorphism_from_B6(G):
	TSS_in_G_size_3 = up_to_conjugation(G,TSS_size_n(G,3))

	possible_homomorphism = false

	for TSS in TSS_in_G_size_3:
		## once again, we can assume without loss of generality that 
		## s1 is sent to TSS[0], s3 to TSS[1] and s5 to TSS[2]

		canidates_s2 = []

		## it must be 1. in the same conjugacy class of s1, 2. not be 
		## in TSS  and 3. braid with s1 and s3, and commute 
		## with s4 and s5

		for elt in G.conjugacy_class(TSS[0]).list():
			if elt in TSS:
				continue
			if braids_with(elt,TSS[0]) and braids_with(elt, TSS[1]) and commutes_with(elt,TSS[2]):
				canidates_s2.append(elt)
				return true
	return false


def print_homomorphism_from_B6(G, name):
	possible_homomorphism = possible_homomorphism_from_B6(G)

	if not possible_homomorphism:
		print("\nThere is no non-cyclic homomorphism B_6 -> " + name + "\n")
	else:
		print("\nThere might be a non-cyclic homomorphism B_6 -> " + name + "\n")


########################################################### 
##                    Maps from B_7'
###########################################################

## Similar to the first function, this can rule out the possibility
## of a homomorphism B_7' -> G

def possible_homomorphism_from_com_B7(G):

	## again, we will use the Lin presentation of B_7'

	## from TSS theory, we know both {c1,c3} and {c2,c4}
	## must be sent to TSS of size 2. 

	TSS_in_G_size_2 = TSS_size_2(G)   
	
	reverse_TSS_in_G = [] ## contents of TSS_in_G_size_2, but reversed
	for TSS in TSS_in_G_size_2:
		reverse_TSS_in_G.append([TSS[1],TSS[0]])

	TSS_including_reversed = TSS_in_G_size_2 + reverse_TSS_in_G

	#print(TSS_in_G_size_2 + reverse_TSS_in_G)

	for TSS1 in up_to_conjugation(G,TSS_in_G_size_2):

		## just as before, we can assume c1,c3 are sent to
		## TSS1[0] and TSS1[1] respectively (by conjugating the map appropriately)
		## however, we cannot make the same assumption of c2 and c4 simultaneously
		## so we also have too loop through the "reverse" of each totally
		## symmetric list in TSS_in_G_size_2. Hence reverse_TSS_in_G.

		c1 = TSS1[0]
		c3 = TSS1[1]

		for TSS2 in TSS_including_reversed:
			c2 = TSS2[0]
			c4 = TSS2[1]

			## the c_i's must braid, commute, and not commute appropriately

			if not braids_with(c1,c2) or not braids_with(c2,c3) or not braids_with(c3,c4):
				continue
			if commutes_with(c1,c2) or commutes_with(c2,c3) or commutes_with(c3,c4):
				continue
			if not commutes_with(c1,c3) or not commutes_with(c1,c4) or not commutes_with(c2,c4):
				continue
			for u in G.conjugacy_class(c2 * c1.inverse()).list():

				## the following are from the Lin relations
				w = u * c1 * u.inverse()     ## 1st relation
				v  = c2.inverse() * u * c2   ## 5th relation (a)

				r2 = (u*w*u.inverse() == w*w*c1.inverse()*w)
				r3 = (v*c1*v.inverse() == c1.inverse()*w)
				r4 = (v*w*v.inverse() == c1.inverse()*w*c1.inverse()*w*c1.inverse()*w*c1.inverse()*c1.inverse()*w)
				r5b = (u*c3 == c3*v)
				r5c = (u*c4 == c4*v)
				r6a = (v*c2 == c2*u.inverse()*v)
				r6b = (v*c3 == c3*u.inverse()*v)
				r6c = (v*c4 == c4*u.inverse()*v)

				if r2 and r3 and r4 and r5b and r6a and r6b:
					return True

			
	return False


def print_homomorphism_from_com_B7(G, name):
	possible_homomorphism = possible_homomorphism_from_com_B7(G)

	if not possible_homomorphism:
		print("\nThere is no non-trivial homomorphism B_7' -> " + name + "\n")
	else:
		print("\nThere might be a non-trivial homomorphism B_7' -> " + name + "\n")



########################################################### 
##                    Maps from B_8'
###########################################################


def possible_homomorphism_from_com_B8(G):

	## again, we will use the Lin presentation of B_8'

	## from TSS theory, we know both {c1,c3} and {c2,c4}
	## must be sent to TSS of size 2. 

	#TSS_in_G_size_2 = TSS_size_2(G)  
	TSS_in_G_size_2 = alternate_TSS_size_2(G) 
	
	reverse_TSS_in_G = [] ## contents of TSS_in_G_size_2, but reversed
	for TSS in TSS_in_G_size_2:
		reverse_TSS_in_G.append([TSS[1],TSS[0]])

	TSS_including_reversed = TSS_in_G_size_2 + reverse_TSS_in_G

	#print(TSS_in_G_size_2 + reverse_TSS_in_G)

	for TSS1 in up_to_conjugation(G,TSS_in_G_size_2):

		## just as before, we can assume c1,c3 are sent to
		## TSS1[0] and TSS1[1] respectively (by conjugating the map appropriately)
		## however, we cannot make the same assumption of c2 and c4 simultaneously
		## so we also have too loop through the "reverse" of each totally
		## symmetric list in TSS_in_G_size_2. Hence reverse_TSS_in_G.

		c1 = TSS1[0]
		c3 = TSS1[1]

		for TSS2 in TSS_including_reversed:
			c2 = TSS2[0]
			c4 = TSS2[1]

			## the c_i's must braid, commute, and not commute appropriately

			if not braids_with(c1,c2) or not braids_with(c2,c3) or not braids_with(c3,c4):
				continue
			if commutes_with(c1,c2) or commutes_with(c2,c3) or commutes_with(c3,c4):
				continue
			if not commutes_with(c1,c3) or not commutes_with(c1,c4) or not commutes_with(c2,c4):
				continue
			for u in G.conjugacy_class(c2 * c1.inverse()).list():

				## the following are from the Lin relations
				w = u * c1 * u.inverse()     ## 1st relation
				v  = c2.inverse() * u * c2   ## 5th relation (a)

				r2 = (u*w*u.inverse() == w*w*c1.inverse()*w)
				r3 = (v*c1*v.inverse() == c1.inverse()*w)
				r4 = (v*w*v.inverse() == c1.inverse()*w*c1.inverse()*w*c1.inverse()*w*c1.inverse()*c1.inverse()*w)
				r5b = (u*c3 == c3*v)
				r5c = (u*c4 == c4*v)
				r6a = (v*c2 == c2*u.inverse()*v)
				r6b = (v*c3 == c3*u.inverse()*v)
				r6c = (v*c4 == c4*u.inverse()*v)

				if r2 and r3 and r4 and r5b and r6a and r6b:
					for c5 in G.conjugacy_class(c1).list():
						if not commutes_with(c5,c1) or not commutes_with(c5,c2) or not commutes_with(c5,c3):
							continue
						if not braids_with(c4,c5):
							continue

						r5d = (u*c5 == c5*v)
						r6d = (v*c5 == c5*u.inverse()*v)

						if r5d and r6d:
							return True

			
	return False


def print_homomorphism_from_com_B8(G, name):
	possible_homomorphism = possible_homomorphism_from_com_B8(G)

	if not possible_homomorphism:
		print("\nThere is no non-trivial homomorphism B_8' -> " + name + "\n")
	else:
		print("\nThere might be a non-trivial homomorphism B_8' -> " + name + "\n")


