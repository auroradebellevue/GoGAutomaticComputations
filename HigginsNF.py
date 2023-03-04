#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 26 13:59:53 2020

@author: aurora 

This file completes all the FSA combining steps to produce the Higgins NF of
a fundamental group of a graph of groups. The FSA's should all have identical
alphabets otherwise these functions will fail.
"""
##############################################################################
# imported files
import os
import subprocess
import FSA_manipulations as FSA

##############################################################################
# global variables
#kbmag_ftn_dir="/home/aurora/gap-4.11.1/pkg/kbmag-1.5.9/standalone/bin/"
#############################################################################
def inflated_higgins_nf(gog, file_dir, kbmag_ftn_dir):
    """
    Parameters
    ----------
    gog : graph of groups

    Returns
    -------
    FSA which accepts the normal forms of the fundamental group of graph of groups

    """ 
    #vertex group word acceptor and edge group coset rep word acceptors 
    for v in gog.vertices:
        v.group.wa=FSA.create_fsa_from_file(file_dir + v.group.wa_file)
        v.group.wa.flags.append("wa")

    for e in gog.edges:
        e.fg.coswa=FSA.create_fsa_from_file(file_dir + e.fg.coswa_file)
        e.rg.coswa=FSA.create_fsa_from_file(file_dir + e.rg.coswa_file)
        e.fg.coswa.flags.append("coswa")
        e.rg.coswa.flags.append("coswa")
    #make a full list of alphabet elements from all vertex groups
    #TODO: test for 2 or more vertices; 2 is tested, 
    full_alpha=[]
    for v in gog.vertices:
        for l in v.group.wa.alphabet.names:
            if l not in full_alpha:
                full_alpha.append(l)
    for e in gog.edges:
        for l in e.fg.coswa.alphabet.names:
            if l not in full_alpha:
                full_alpha.append(l)
        for l in e.rg.coswa.alphabet.names:
            if l not in full_alpha:
                full_alpha.append(l)
    #add letters that are missing when comparing to full_alpha and add empty transitions for the new letters
    #calculate the permutation to organize each edge and vertex alphabet
    #to match the order in full_alpha. The permutation will be written using 
    #the following rule: the number in position i indicates the 
    #new index of the ith letter of the alphabet, etc.
    #Then call permute_order
    for v in gog.vertices:
        v_alphabet=v.group.wa.alphabet.names
        for l in full_alpha:
            if l not in v_alphabet:
                v.group.wa.new_letter(l)
        perm=[]
        for l in v_alphabet:
             perm.append(full_alpha.index(l))
        FSA.permute_order(v.group.wa, perm)
    
    for e in gog.edges:
        ef_alphabet=e.fg.coswa.alphabet.names
        for l in full_alpha:
            if l not in ef_alphabet:
                e.fg.coswa.new_letter(l)
        perm=[]
        for l in ef_alphabet:
             perm.append(full_alpha.index(l))
        FSA.permute_order(e.fg.coswa, perm)
        er_alphabet=e.rg.coswa.alphabet.names
        for l in full_alpha:
            if l not in er_alphabet:
                e.rg.coswa.new_letter(l)
        perm=[]
        for l in er_alphabet:
             perm.append(full_alpha.index(l))
        FSA.permute_order(e.rg.coswa, perm)
        
    #add stable letter and inverse to the end of the alphabet of each fsa 
    stable_letter_list=[]
    for e in gog.edges: 
        stable_letter_list.append(e.stab_ltr[0])
    for e in gog.edges:
        for t in stable_letter_list:
            e.fg.coswa.new_stab_letter(t)
            e.rg.coswa.new_stab_letter(t)
        e.fg.coswa.print_fsa(file_dir + "full_alpha_"+e.fg.coswa_file)
        e.rg.coswa.print_fsa(file_dir + "full_alpha_"+e.rg.coswa_file)
    for v in gog.vertices:
        for t in stable_letter_list:
            v.group.wa.new_stab_letter(t)
        v.group.wa.print_fsa(file_dir + "full_alpha_"+ v.group.wa_file)
    
    print("Attempted to add extra letters to fsa's and print back to file")
    #build coset blocks for each edge (and reverse edge) 
    for e in gog.edges:
        e.coset_block_files=build_coset_blocks(e, file_dir, kbmag_ftn_dir)
    
    #build product of coset blocks for each edge 
    or_coset_block_files=[]
    for e in gog.edges:
        subprocess.run([kbmag_ftn_dir+"fsaor",
                        file_dir + e.coset_block_files[0],
                        file_dir + e.coset_block_files[1],
                        file_dir + "OR_coset_blocks_e" + e.label])
        or_coset_block_files.append("OR_coset_blocks_e"+e.label)
    #Build Language accepting the words in any one of the coset blocks (take "or" of all coset blocks) 
    if len(or_coset_block_files)>=2:
        subprocess.run([kbmag_ftn_dir+"fsaor",
                        file_dir + or_coset_block_files[0],
                        file_dir + or_coset_block_files[1],
                        file_dir + "ALL_coset_blocks"])
        if len(or_coset_block_files)>2:
            for block_file in or_coset_block_files[2:]:
                subprocess.run([kbmag_ftn_dir+"fsaor",
                        file_dir + block_file,
                        file_dir + "ALL_coset_blocks",
                        file_dir + "ALL_coset_blocks"])
    elif len(or_coset_block_files)==1:
        os.rename(file_dir + "OR_coset_blocks_e0", 
                  file_dir + "ALL_coset_blocks")
    else:
        print("something went wrong making the coset blocks since there are no edges?")
    #Take the star of the OR language of all coset blocks #step is GOG ready
    subprocess.run([kbmag_ftn_dir+"fsastar", file_dir+ "ALL_coset_blocks"]) 
    print("coset blocks are made")
    #Inflated language of normal forms with "bad words"
    subprocess.run([kbmag_ftn_dir+"fsaconcat", 
                    file_dir + "full_alpha_"+ gog.vertices[0].group.wa_file,
                    file_dir + "ALL_coset_blocks.star",
                    file_dir + "infl_nf_with_bad_words"])
    print("inflated normal forms with bad words fsa is printed to file")
    bad_word_collection_files=[]
    #Collection 1 of bad words: Language with the subwords stable letter followed by its inverse 
    #alphabet.wa.star is created in file_dir in the process, and can be accessed in other steps
    print("collection 1 of bad words started")
    tT_subword_files=[]
    bad_word_collection_files.append("ALL_Coll1")
    for e in gog.edges:
        e.lL_subword_file=ltr_Ltr_subword(e, file_dir, kbmag_ftn_dir)
        tT_subword_files.append(e.lL_subword_file)
    if len(tT_subword_files)>=2:
        subprocess.run([kbmag_ftn_dir+"fsaor",
                        file_dir + tT_subword_files[0],
                        file_dir + tT_subword_files[1],
                        file_dir + "ALL_Coll1"])
        if len(tT_subword_files)>2:
            for tT_file in tT_subword_files[2:]:
                subprocess.run([kbmag_ftn_dir+"fsaor",
                        file_dir + tT_file,
                        file_dir + "ALL_Coll1",
                        file_dir + "ALL_Coll1"])
    elif len(tT_subword_files)==1:
        os.rename(file_dir + tT_subword_files[0], file_dir + "ALL_Coll1")
    else:
        print("Collection 1 of bad words not needed or something went wrong")
        bad_word_collection_files[0]="None"
    #Collection 2 of bad words: sequence of stable letters which do not form a path on the graph with some 
    #valid coset representive in between the stable letters
    print("Collection 2 of bad words started")
    total=len(gog.edges)
    invalid_next_stable_letter=[] #entry i is a list for the invalid followers of t_i then the invalid followers of T_i
    for first_e in gog.edges: #TODO this part depends on the edges being listed with index in bijection with their label 
        temp_list=[[],[], first_e.label]
        #print(first_e.label)
        for next_e in gog.edges:
            if next_e.iv != first_e.tv:
                temp_list[0].append("t"+ next_e.label)
            if next_e.tv != first_e.tv:
                temp_list[0].append("T"+ next_e.label)
            if next_e.iv != first_e.iv: 
                temp_list[1].append("t"+ next_e.label)
            if next_e.tv != first_e.iv:
                temp_list[1].append("T"+ next_e.label)
        invalid_next_stable_letter.append(temp_list)
    #print(invalid_next_stable_letter)
    invalid_stable_letter_pair_files=[] #list all files for coset blocks followed by an invalid stable letter
    for e in gog.edges:
        alpha=e.fg.coswa.alphabet
        for entry in invalid_next_stable_letter:
            if e.label in entry:
                index1=invalid_next_stable_letter.index(entry)
        #print(e.label, index1)
        row=invalid_next_stable_letter[index1]
        for letter in row[0]:
            next_e_fsa=FSA.create_small_fsa(alpha, letter)
            next_e_fsa.print_fsa(file_dir+letter)
            filename=str("badpath_t"+e.label+"_"+letter)
            subprocess.run([kbmag_ftn_dir+"fsaconcat",
                file_dir + e.coset_block_files[0],
                file_dir + letter,
                file_dir + filename])
            invalid_stable_letter_pair_files.append(filename)
        for letter in row[1]:
            next_e_fsa=FSA.create_small_fsa(alpha, letter)
            next_e_fsa.print_fsa(file_dir+letter)
            filename=str("badpath_T"+e.label+"_"+letter)
            subprocess.run([kbmag_ftn_dir+"fsaconcat",
                file_dir + e.coset_block_files[1],
                file_dir + letter,
                file_dir + filename])
            invalid_stable_letter_pair_files.append(filename)
    #build the language which has a bad pair as a subword 
    subwords_invalid_stable_letter_pair_files=[]
    bad_word_collection_files.append("ALL_Coll2")
    for x in invalid_stable_letter_pair_files:
        subprocess.run([kbmag_ftn_dir+"fsaconcat",
                file_dir + "alphabet.wa.star",
                file_dir + x,
                file_dir + "AW_" + x]) #AW is short for anyword. kbmag gets mad when filenames get long
        subprocess.run([kbmag_ftn_dir+"fsaconcat",
                file_dir + "AW_" + x,
                file_dir + "alphabet.wa.star",
                file_dir + "AW_" + x + "_AW"])
        subwords_invalid_stable_letter_pair_files.append("AW_"+ x + "_AW")
    #Take "or" operation of all of them
    if len(subwords_invalid_stable_letter_pair_files)>=2:
        subprocess.run([kbmag_ftn_dir+"fsaor",
                        file_dir + subwords_invalid_stable_letter_pair_files[0],
                        file_dir + subwords_invalid_stable_letter_pair_files[1],
                        file_dir + "ALL_Coll2"])
        if len(subwords_invalid_stable_letter_pair_files)>2:
            for file in subwords_invalid_stable_letter_pair_files[2:]:
                subprocess.run([kbmag_ftn_dir+"fsaor",
                                file_dir + file,
                                file_dir + "ALL_Coll2",
                                file_dir + "ALL_Coll2"])
    elif len(subwords_invalid_stable_letter_pair_files)==1:
        os.rename(file_dir + subwords_invalid_stable_letter_pair_files[0], 
                  file_dir + "ALL_Coll2")
    else:
        print("Collection 2 of bad words not needed or something went wrong")
        bad_word_collection_files[1]="None"
    #Collection 3 of bad words: initial stable letter does not correspond to an 
    #edge whose initial vertex is v0
    #find all vertices which don't have v0 as the initial (then t cannot follow)
    #find all vertices which don't have v0 as the terminal (then T cannot follow)
    print("Collection 3 of bad words started")
    invalid_first_stable_letter=[]
    for e in gog.edges:
        if e.iv != 0:
            invalid_first_stable_letter.append("t"+e.label)
        if e.tv != 0:
            invalid_first_stable_letter.append("T"+e.label)
    #build language (v0 normal forms)t(anyword)
    invalid_first_stable_letter_subword_files=[]
    for t in invalid_first_stable_letter:
        invalid_t_fsa=FSA.create_small_fsa(alpha, t)
        invalid_t_fsa.print_fsa(file_dir+t+".wa")
        subprocess.run([kbmag_ftn_dir+"fsaconcat",
                file_dir + "full_alpha_v0gp.wa",
                file_dir + t + ".wa",
                file_dir + "v0_NF_" + t + ".wa"]) #NF is normal forms   
        subprocess.run([kbmag_ftn_dir+"fsaconcat",
                file_dir + "v0_NF_" + t + ".wa",
                file_dir + "alphabet.wa.star",
                file_dir + "v0_NF_" + t + "_AW" + ".wa"]) #AW is anyword
        invalid_first_stable_letter_subword_files.append("v0_NF_" + t + "_AW" + ".wa")
    #Take or statement of the languages whose first stable letter is invalid
    bad_word_collection_files.append("ALL_Coll3")
    if len(invalid_first_stable_letter_subword_files)>=2:
        subprocess.run([kbmag_ftn_dir+"fsaor",
                        file_dir + invalid_first_stable_letter_subword_files[0],
                        file_dir + invalid_first_stable_letter_subword_files[1],
                        file_dir + "ALL_Coll3"])
        if len(invalid_first_stable_letter_subword_files)>2:
            for file in invalid_first_stable_letter_subword_files[2:]:
                subprocess.run([kbmag_ftn_dir+"fsaor",
                        file_dir + file,
                        file_dir + "ALL_Coll3",
                        file_dir + "ALL_Coll3"])
    elif len(invalid_first_stable_letter_subword_files)==1:
        os.rename(file_dir + invalid_first_stable_letter_subword_files[0], 
                  file_dir + "ALL_Coll3")
    else:
        print("Collection 3 of bad words not needed or something went wrong")
        bad_word_collection_files[2]="None"
    
    #Collection 4: build the language that ends with a stable letter that 
    #corresponds to any direction of an edge on the spanning tree
    #use spanning tree from Graph of Groups
    print("Collection 4 of bad words started")
    tree=gog.spanning_tree
    if tree !=[]:
        tree_stable_letters=[]
        for e in gog.edges:
            if "e"+e.label in tree:
                tree_stable_letters.extend(["t"+e.label, "T"+e.label])
            if "e"+e.label+"r" in tree:
                tree_stable_letters.extend(["t"+e.label, "T"+e.label])
        #make fsa that accepts a letter from the list of tree stable letters
        alpha=gog.vertices[0].group.wa.alphabet
        sts=FSA.states(2, [1], [2])
        row1=[]
        row2=[]
        for a in alpha.names:
            if a in tree_stable_letters:
                row1.append(2)
                row2.append(0)
            else: 
                row1.append(0)
                row2.append(0)
            tr=[row1, row2]
            t=FSA.table("dense deterministic", tr)
            tree_fsa=FSA.fsa(alpha, sts, t)
            tree_fsa.print_fsa(file_dir + "tree_stable_letters")
        #concatenate the language of any word before the tree stable letters
        subprocess.run([kbmag_ftn_dir+"fsaconcat",
                        file_dir + "alphabet.wa.star",
                        file_dir + "tree_stable_letters",
                        file_dir + "ALL_Coll4"])
        bad_word_collection_files.append("ALL_Coll4")
    else:
        print("Collection 4 of bad words not needed")
        bad_word_collection_files.append("None")
    #Or of all collections:
    if bad_word_collection_files[0]!="None":
        os.system('cp '+ file_dir+bad_word_collection_files[0]+" "+ file_dir+'ALL_Colls')
        #os.rename(file_dir + bad_word_collection_files[0], 
        #          file_dir + "ALL_Colls")
    if bad_word_collection_files[1]!="None":
        subprocess.run([kbmag_ftn_dir+"fsaor",
                file_dir + "ALL_Colls",
                file_dir + bad_word_collection_files[1],
                file_dir + "ALL_Colls"])
    if bad_word_collection_files[2]!="None":
        subprocess.run([kbmag_ftn_dir+"fsaor",
                file_dir + "ALL_Colls",
                file_dir + bad_word_collection_files[2],
                file_dir + "ALL_Colls"])
    if bad_word_collection_files[3]!="None":
        subprocess.run([kbmag_ftn_dir+"fsaor",
                file_dir + "ALL_Colls",
                file_dir + bad_word_collection_files[3],
                file_dir + "ALL_Colls"])
    #Inflated language of normal forms without "bad words"
    subprocess.run([kbmag_ftn_dir+"fsaandnot",
                file_dir + "infl_nf_with_bad_words",
                file_dir + "ALL_Colls",
                file_dir + "InflatedHigginsNF.wa"])
    return FSA.create_fsa_from_file(file_dir + "InflatedHigginsNF.wa")
    
def ltr_Ltr_subword(e, file_dir, kbmag_ftn_dir):
    """
    Parameters
    ----------
    e : edge object
    file_dir : TYPE
        DESCRIPTION.

    Returns
    -------
    None. Builds FSA that accepts the language of words which have a subword
    of the form stable letter adjacent to its inverse. 

    """
    #make smaller fsa's 
    stab_wa=FSA.create_small_fsa(e.fg.coswa.alphabet, e.stab_ltr[0])
    stab_wa.print_fsa(file_dir + "stab" + e.label+ ".wa") 
    stab_inv_wa=FSA.create_small_fsa(e.rg.coswa.alphabet, e.stab_ltr[1])
    stab_inv_wa.print_fsa(file_dir + "stab_inv" + e.label + ".wa")
    
    subprocess.run([kbmag_ftn_dir+"fsaconcat", 
                    file_dir + "stab" + e.label +".wa",
                    file_dir + "stab_inv" + e.label + ".wa",
                    file_dir + e.stab_ltr[0] + e.stab_ltr[1] + ".wa"])
    subprocess.run([kbmag_ftn_dir+"fsaconcat", 
                    file_dir + "stab_inv" + e.label + ".wa",
                    file_dir + "stab"+ e.label + ".wa",
                    file_dir + e.stab_ltr[1] + e.stab_ltr[0] + ".wa"])
    #inflated alphabet kleene star
    e.fg.coswa.alphabet.create_star_fsa("alphabet.wa", file_dir, kbmag_ftn_dir) #this must be called after the alphabet has been inflated
    
    #concatenate these to create the "lL" subwords
    subprocess.run([kbmag_ftn_dir+"fsaconcat",
                file_dir + "alphabet.wa.star",
                file_dir + e.stab_ltr[0] + e.stab_ltr[1] + ".wa",
                file_dir + "anyword_" + e.stab_ltr[0] + e.stab_ltr[1] + ".wa"])
    subprocess.run([kbmag_ftn_dir+"fsaconcat",
                file_dir + "anyword_" + e.stab_ltr[0] + e.stab_ltr[1] + ".wa",
                file_dir + "alphabet.wa.star",
                file_dir + "anyword_" + e.stab_ltr[0] + e.stab_ltr[1] + "_anyword.wa"])
    subprocess.run([kbmag_ftn_dir+"fsaconcat",
                file_dir + "alphabet.wa.star",
                file_dir + e.stab_ltr[1] + e.stab_ltr[0] + ".wa",
                file_dir + "anyword_" + e.stab_ltr[1] + e.stab_ltr[0] + ".wa"])
    subprocess.run([kbmag_ftn_dir+"fsaconcat",
                file_dir + "anyword_" + e.stab_ltr[1] + e.stab_ltr[0] + ".wa",
                file_dir + "alphabet.wa.star",
                file_dir + "anyword_" + e.stab_ltr[1] + e.stab_ltr[0] + "_anyword.wa"])
    subprocess.run([kbmag_ftn_dir+"fsaor",
                file_dir + "anyword_" + e.stab_ltr[0] + e.stab_ltr[1] + "_anyword.wa",
                file_dir + "anyword_" + e.stab_ltr[1] + e.stab_ltr[0] + "_anyword.wa",
                file_dir + "stable_lL_e" + e.label])
    return "stable_lL_e" + e.label      

def build_coset_blocks(e, file_dir, kbmag_ftn_dir):
    """
    Parameters
    ----------
    e : edge object

    Returns
    -------
    List of names of new FSA files. 
    FSAs that accept the languages t(right cosets of G_e) and 
    T(right cosets of G_(e reverse)) are saved as files.
    """
    #prepare files in the folder

    print("In build_coset blocks function")
    stab_wa=FSA.create_small_fsa(e.fg.coswa.alphabet, e.stab_ltr[0])
    stab_wa.print_fsa(file_dir + "stab" + str(e.label) +".wa") 
    stab_inv_wa=FSA.create_small_fsa(e.rg.coswa.alphabet, e.stab_ltr[1])
    stab_inv_wa.print_fsa(file_dir + "stab_inv" + str(e.label) +".wa")
    
    #concatenate the stable letter to the beginning of the coset language
    subprocess.run([kbmag_ftn_dir+"fsaconcat",
                    file_dir+"stab" + str(e.label) + ".wa",
                    file_dir+"full_alpha_" + e.fg.coswa_file,
                    file_dir+"coset_block_e"+ str(e.label) + "-f"])
    subprocess.run([kbmag_ftn_dir+"fsaconcat",
                    file_dir+"stab_inv"+str(e.label)+".wa",
                    file_dir+"full_alpha_" + e.rg.coswa_file,
                    file_dir+"coset_block_e"+ str(e.label) + "-r"])

    return ["coset_block_e"+ str(e.label) + "-f", "coset_block_e"+ str(e.label) + "-r"]

def deflation(gog, file_dir):
    """From the fsa accepting the inflated language of normal forms, calculate
    the language of normal forms by taking the deflation with respect to the
    tree. The deflation map sends each stable letter corresonding to an 
    edge on the tree to the empty word and leaves the rest of the alphabet alone."""
    tree=gog.spanning_tree    
    #print(tree)
    tree_stable_letters=[]
    for e in gog.edges:
        #print(e.label)
        #print("e"+e.label, "e"+e.label+'r')
        if "e"+e.label in tree:
            tree_stable_letters.append("t"+e.label)
            tree_stable_letters.append("T"+e.label)
        if "e"+e.label+"r" in tree:
            tree_stable_letters.append("t"+e.label)
            tree_stable_letters.append("T"+e.label)
    temp_fsa=gog.inflated_wa
    if tree_stable_letters==[]:
        print("No deflation required. gog.wa and gog.inflated_wa store the FSA accepting the normal forms")
    else:
        #make all tree letters the empty word
        print(tree_stable_letters)
        for t in tree_stable_letters:
            index_t=temp_fsa.alphabet.names.index(t)
            temp_fsa.alphabet.names[index_t]="_"
            temp_fsa.alphabet.size+=-1
        temp_fsa.alphabet.size+=1 #to adjust for '_'
        #make inflated word acceptor non-det fsa since there are epsilon transitions now
        temp_fsa=FSA.make_ndfsa(temp_fsa)
        print(temp_fsa.states.size, temp_fsa.states.accepting)
        #print(temp_fsa.table.transitions)
        #make temp_fsa which is nd with epsilon into e free ndfsa
        temp_fsa=FSA.e_free_ndfsa(temp_fsa) 
        print(temp_fsa.states.size, temp_fsa.states.accepting)
        #print(temp_fsa.table.transitions)
        #change temp_fsa from e free ndfsa to dfsa
        temp_fsa=FSA.make_det(temp_fsa)
        print(temp_fsa.states.size, temp_fsa.states.accepting)
        print("Deflation of inflated normal forms complete. The FSA accepting the normal forms is stored in gog.wa")
    
    return temp_fsa
        #TODO test and output