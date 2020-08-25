
## This file is ONLY FOR COMPUTATIONS



########################################################### 
##                   Sage Import Shenanigans
###########################################################

# import os
# os.system('sage --preparse TotallySymmetricSets.sage')
# os.system('mv TotallySymmetricSets.sage.py TotallySymmetricSets.py')
# from TotallySymmetricSets import *

import os
os.system('sage --preparse HomomorphismCheckingAlgorithms.sage')
os.system('mv HomomorphismCheckingAlgorithms.sage.py HomomorphismCheckingAlgorithms.py')
from HomomorphismCheckingAlgorithms import *

import time



########################################################### 
##                         n=6 
###########################################################

print("\n------------  n=6  ------------")
print("Expected time: 10 seconds")

start = time.time()

print_homomorphism_from_com_B6(PSL27, "PSL(2,7)")

print_homomorphism_from_B6(PGL29, "PGL(2,9)")

print_homomorphism_from_B6(M_10, "Mathieu 10")

end = time.time()

print(end-start)


########################################################### 
##                         n=7 
###########################################################


print("\n------------  n=7  ------------")
print("Expected time: 2.5 minutes")

start = time.time()

print_homomorphism_from_com_B7(PSL28, "PSL(2,8)")

print_homomorphism_from_com_B7(PSL211, "PSL(2,11)")

print_homomorphism_from_com_B7(PSL213, "PSL(2,13)")

print_homomorphism_from_com_B7(PSL217, "PSL(2,17)")

print_homomorphism_from_com_B7(PSL216, "PSL(2,16)")

end = time.time()

print(end-start)

########################################################### 
##                         n=8 
###########################################################


print("\n------------  n=8  ------------")
print("Expected time: 50 minutes")

start = time.time()

print_homomorphism_from_com_B8(PSL33,"PSL(3,3)")

print_homomorphism_from_com_B8(G22_prime, "G_2(2)'")

print_homomorphism_from_com_B8(M_11,"Mathieu 11")

print_homomorphism_from_com_B8(PSL34, "PSL(3,4)")

end = time.time()

print(end-start)


# ########################################################### 
# ##                       Testing 
# ###########################################################

# print("\n------------  Testing  ------------")
# print("Expected time: 40 minutes")

# start = time.time()

# print_homomorphism_from_com_B6(A6, "A_6")

# print_homomorphism_from_B6(S6, "S_6")

# print_homomorphism_from_com_B7(A7, "A_7")

# print_homomorphism_from_com_B8(A8, "A_8")


# end = time.time()

# print(end-start)