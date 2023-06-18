#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 17 19:21:29 2023

@author: aurora
"""
import FSA_manipulations as FSA
import subprocess
verbose = True
file_dir = "/home/aurora/ArtinTest/"
kbmag_ftn_dir = "/usr/share/gap/pkg/kbmag-1.5.9/standalone/bin/"

artinWAFile = "/home/aurora/ArtinEx2/OutputFolder/Pi1GraphofGroupsWordAcceptor.wa"


def testMultiplier(testFSA, fsaString):
    eps_alpha = testFSA.alphabet.base.copy()
    eps_alpha.append("_")
    letter_alphabet = FSA.alphabet(len(eps_alpha), eps_alpha, "identifiers")
    #find initial state
    if fsaString == 'IdWord':
        init_state = {1}
    else:
        init_state = {}
    #The intial state is always an accept state since we are not checking the 
    #pairs (_, x) for x in X
    test_acc_states = [{1}]

    test_table=[]
    isAcceptState = True
    #find the transition table for the projection onto the first coordinate
    for iRow in range(len(testFSA.table.transitions)):
        temp_row = []
        for iFirstLetter in range(len(testFSA.alphabet.base)+1):
            new_st = set()
            for iSecondLetter in range(len(testFSA.alphabet.base)+1):
                alpha_index = iFirstLetter*(len(testFSA.alphabet.base)+1)+iSecondLetter
                if alpha_index != len(eps_alpha)*len(eps_alpha)-1:
                    new_st.add(testFSA.table.transitions[iRow][alpha_index])
            new_st = FSA.rem_0(new_st)
            #check if it's an accept state
            for iState in new_st:
                if iState not in testFSA.states.accepting:
                    isAcceptState = False
            if isAcceptState == True:
                test_acc_states.append(new_st)
            #reset to True
            isAcceptState = True
            temp_row.append(new_st)
        test_table.append(temp_row)


    #print("test table info", len(test_table))
    #print("new test table", test_table)   
    #print("accept states", test_acc_states)
    #define state object
    letter_states = FSA.states(len(test_table), [{1}], test_acc_states)
    #define table object
    letter_table_object = FSA.table("non-det dense", test_table)
    #define the ndfsa with epsilon transitions
    letter_ndfsa_with_eps = FSA.ndfsa(letter_alphabet, letter_states, letter_table_object)
    #make it an ndfsa without epsilon transitions
    letter_ndfsa = FSA.e_free_ndfsa(letter_ndfsa_with_eps)
    #make it a det fsa 
    letter_fsa = FSA.make_det(letter_ndfsa)
    #print("final table", letter_fsa.table.transitions)
    #print("final state list", letter_fsa.states.size)
    #print("final accept states", letter_fsa.states.accepting)
    #print out fsa file and compare the projection to the NF
    projection_file = "proj1-Artin2-Letter-"+fsaString+".txt"
    letter_fsa.print_fsa(file_dir + projection_file)
    output = subprocess.run([kbmag_ftn_dir+"fsalequal", artinWAFile, file_dir+projection_file])
    print("fsaequal output for letter: "+fsaString+" "+ str(output.returncode)+ "\n")

#Test k=2
print("Testing k=2 multipliers")
artin2_id = FSA.create_fsa_from_file("/home/aurora/ArtinEx2/OutputFolder/fsaminLetter_IdWord.txt")
artin2_a = FSA.create_fsa_from_file("/home/aurora/ArtinEx2/OutputFolder/fsaminLetter_a.txt")
artin2_A = FSA.create_fsa_from_file("/home/aurora/ArtinEx2/OutputFolder/fsaminLetter_A.txt")
artin2_b = FSA.create_fsa_from_file("/home/aurora/ArtinEx2/OutputFolder/fsaminLetter_b.txt")
artin2_B = FSA.create_fsa_from_file("/home/aurora/ArtinEx2/OutputFolder/fsaminLetter_B.txt")
artin2_c = FSA.create_fsa_from_file("/home/aurora/ArtinEx2/OutputFolder/fsaminLetter_c.txt")
artin2_C = FSA.create_fsa_from_file("/home/aurora/ArtinEx2/OutputFolder/fsaminLetter_C.txt")
artin2_t = FSA.create_fsa_from_file("/home/aurora/ArtinEx2/OutputFolder/fsaminLetter_t0.txt")
artin2_T = FSA.create_fsa_from_file("/home/aurora/ArtinEx2/OutputFolder/fsaminLetter_T0.txt")
#Test the IdWord
testMultiplier(artin2_id, "IdWord")
#Test the a
testMultiplier(artin2_a, "a")
#Test the A
testMultiplier(artin2_A, "A")
#Test the b
testMultiplier(artin2_b, "b")
#Test the B
testMultiplier(artin2_B, "B")
#Test the c
testMultiplier(artin2_c, "c")
#Test the C
testMultiplier(artin2_C, "C")
#Test the t0
testMultiplier(artin2_t, "t0")
#Test the T0
testMultiplier(artin2_T, "T0")

#Test k=3
print("Testing k=3 multipliers")
artin3_id = FSA.create_fsa_from_file("/home/aurora/ArtinEx3/OutputFolder/fsaminLetter_IdWord.txt")
artin3_a = FSA.create_fsa_from_file("/home/aurora/ArtinEx3/OutputFolder/fsaminLetter_a.txt")
artin3_A = FSA.create_fsa_from_file("/home/aurora/ArtinEx3/OutputFolder/fsaminLetter_A.txt")
artin3_b = FSA.create_fsa_from_file("/home/aurora/ArtinEx3/OutputFolder/fsaminLetter_b.txt")
artin3_B = FSA.create_fsa_from_file("/home/aurora/ArtinEx3/OutputFolder/fsaminLetter_B.txt")
artin3_c = FSA.create_fsa_from_file("/home/aurora/ArtinEx3/OutputFolder/fsaminLetter_c.txt")
artin3_C = FSA.create_fsa_from_file("/home/aurora/ArtinEx3/OutputFolder/fsaminLetter_C.txt")
artin3_t = FSA.create_fsa_from_file("/home/aurora/ArtinEx3/OutputFolder/fsaminLetter_t0.txt")
artin3_T = FSA.create_fsa_from_file("/home/aurora/ArtinEx3/OutputFolder/fsaminLetter_T0.txt")
#Test the IdWord
testMultiplier(artin3_id, "IdWord")
#Test the a
testMultiplier(artin3_a, "a")
#Test the A
testMultiplier(artin3_A, "A")
#Test the b
testMultiplier(artin3_b, "b")
#Test the B
testMultiplier(artin3_B, "B")
#Test the c
testMultiplier(artin3_c, "c")
#Test the C
testMultiplier(artin3_C, "C")
#Test the t0
testMultiplier(artin3_t, "t0")
#Test the T0
testMultiplier(artin3_T, "T0")