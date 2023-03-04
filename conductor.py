#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 16:11:59 2020

@author: Aurora Marks

This python script is the main file to "conduct" and control the processes of 
the files GraphofGroups.py, run-auto.py, HigginsNF.py, and FSA_manipulations.py.

There are four main stages of this procedure: 
    I. Manage the input and internal storage of the graph of groups
    II. Run the appropriate kbmag programs to obtain the automatic structures
    III. Form the Word Acceptor for Normal Forms of the Fundamental Group of
        the graph of groups. 
    IV. Form the Multiplier Automaton
    
    Notes on each stage:
        I. This python script has one input, the absolute path of the directory
        that contains all group rewriting files, subgroup rewriting files, 
        graph.txt which contains the directed adjacency matrix, and 
        file-info.txt which stores the file names. 
        
"""
#import libraries
import sys
import os

#import other scripts
import GraphofGroups as GraphofGroups
import run_auto as runauto
import HigginsNF as HigginsNF
import Multiplier as Multiplier

def unpack_gog(file_info, adj):
    """
    

    Parameters
    ----------
    file_info : TYPE
        DESCRIPTION.
    adj : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    vfiles=[]
    efiles=[]
    isoms=[]
    init_term_vs=[]
    for i in range(len(adj)):
        for j in range(len(adj)):
            #if i=j then the vertex info is in this entry as well as any edges that
            #start and end at the same vertex
            if i==j:
                file_info[i][j]=file_info[i][j].split(',')
                #vertex file with independent ordering
                vfiles.append(file_info[i][j][0].replace(" ", ""))
                for k in range(adj[i][j]):
                    #the subgroup file and the isomorphic subgroup file are saved in efiles
                    #for each directed edge, there are 5 pieces of information
                    #viej-f, viej-f.sub, viej-r, viej-r.sub, isom j
                    efiles.append([file_info[i][j][5*k+1].replace(" ", ""),
                                   file_info[i][j][5*k+2].replace(" ", ""),
                                   file_info[i][j][5*k+3].replace(" ", ""),
                                   file_info[i][j][5*k+4].replace(" ", "")
                                   ])
                    isoms.append(file_info[i][j][5*k+5].replace(" ", ""))
                    init_term_vs.append([i,j])
                    #print(isoms)
            #if i and j differ then there is only subgroup information
            else:
                for k in range(adj[i][j]):
                    file_info[i][j]=file_info[i][j].split(',')
                    #the subgroup file and the isomorphic subgroup file are saved in efiles
                    #for each directed edge, there are 5 pieces of information
                    #viej-f, viej-f.sub, viej-r, viej-r, isom j
                    efiles.append([file_info[i][j][5*k].replace(" ", ""),
                                   file_info[i][j][5*k+1].replace(" ", ""),
                                   file_info[i][j][5*k+2].replace(" ", ""),
                                   file_info[i][j][5*k+3].replace(" ", "")
                                   ])
                    isoms.append(file_info[i][j][5*k+4].replace(" ", ""))
                    init_term_vs.append([i,j])
                    #print(isoms)
    #change isomorphisms to a list of dictionaries
    for i in range(len(isoms)):
        l=int(len(isoms[i])/2)
        x=isoms[i][0:l] + isoms[i][l:l+l]
        y=isoms[i][l:l+l] + isoms[i][0:l]
        table=isoms[i].maketrans(x,y)
        isoms[i]=table
    return vfiles, efiles, isoms, init_term_vs                 

def main(argv):
    """
    Full graph of groups ready until Normal forms
    Parameters
    ----------
    argv : absolute path of the directory that contains the group files

    Returns
    -------
    Something evenutally

    """
    # Upack Input and make Output Folder
    file_dir=argv[0]
    kbmag_ftn_dir=argv[1]
    out_dir = os.path.join(file_dir, r'OutputFolder/')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    # Find files containing graph of groups structure and perform basic processing
    file_input=runauto.find_files(dir_path=file_dir, pattern='file-info.txt', flag='start')
    graph_input=runauto.find_files(dir_path=file_dir, pattern='graph.txt', flag='start')
    travel_input=runauto.find_files(dir_path=file_dir, pattern='fellow-traveling-info.txt', flag='start')
    file_info=runauto.basic_processing(file_input, 'files')
    graph_info=runauto.basic_processing(graph_input, 'graph')
    travel_info=runauto.basic_processing(travel_input, 'travel')
    
    # Convert graph_info into an adjacency matrix with integer entries
    adj, spanning_tree = runauto.get_graph_and_tree(graph_info)
    # Use adjaceny matrix, adj, to unpack and store information from the files
    vfiles, efiles, isoms, init_term_vs=unpack_gog(file_info, adj)
    
    # Testing output:
    print("file info", file_info, type(file_info), len(file_info))
    print("graph info", graph_info, type(graph_info), len(graph_info))
    print(adj)
    print("vertex files" , vfiles)
    print("edge files" , efiles)
    print("isomorphisms stored as ascii code" , isoms)
    print("initial and terminal vertices of corresponding edge", init_term_vs)
    #print (gog)
    
    # Internal Storage
    gog=GraphofGroups.make_gog(vfiles, efiles, init_term_vs, isoms, adj, 
                               file_dir)
    gog.spanning_tree=spanning_tree
    print("finished making Graph of Groups")
    #Create files for other orderings on each vertex
    #USER MUST CREATE EXTRA FILES
    #GraphofGroups.other_orderings(gog)
    
    # Run kbmag programs
    runauto.run_auts(gog,vfiles,efiles,init_term_vs,file_dir,kbmag_ftn_dir)  
    
    
    # Word Acceptor
    gog.inflated_wa=HigginsNF.inflated_higgins_nf(gog, file_dir, kbmag_ftn_dir)
    print("finished making inflated normal forms")
    print("attempting deflation")
    gog.wa=HigginsNF.deflation(gog, file_dir)
    gog.wa.print_fsa(out_dir+gog.wa_file)
    print("Alphabet and transition table of Word Acceptor")
    print(gog.wa.alphabet.names)
    print(gog.wa.table.transitions)
    print("Number of states in word acceptor", gog.wa.states.size)
    
    #Cascade and Free Reduction Testing
    
    """
    test_words = ["bt0baT0bat0ba", 
                  "bt0abT0ba", 
                  "t0abT0bat0", 
                  "t0abT0ba", 
                  "t0aat0bb",
                  "t0T0", 
                  "IdWord",
                  str(),
                  "Bt0T0b", 
                  "t0ABT0t0baT0aAAaBbbBbBbB", 
                  "t0T0bat0a",
                  "T0bbt0a",
                  "bat0a"]
    print('\nTesting cascade and free reduction combination function\n')
    test_result = []
    for i in range(len(test_words)):
        test_result.append(Multiplier.casc_fred_combo(test_words[i], gog, file_dir, kbmag_ftn_dir))
    for i in range(len(test_words)):
        print("Input: ", test_words[i], " Output: ", test_result[i] +"\n")
    """
    
    #Multiplier Automaton
    print("***********************************\n\n")
    print("attempting multiplier building")
    gog.gm=Multiplier.make_gm(gog, travel_info, file_dir, kbmag_ftn_dir, out_dir, False)
    for iLetter in range(len(gog.mult_list)):
        gog.mult_list[iLetter][1].print_fsa(out_dir+"fsaminLetter_"+gog.mult_list[iLetter][0]+".txt")
    print("***********************************\n\n")
    print("The automatic structure is built. The files are called Pi1GraphofGroupsWordAcceptor.wa and fsaminLetter_***.txt")
    print("***********************************\n\n")
if __name__=="__main__":
    main(sys.argv[1:]) #the first argument is the python script

#############################################################################
#testing ground
#file_lines=[["v0, e0-f.sub, e0-r.sub, aAaA","e1-f.sub, e1-r.sub, bBbB ",""],
#["","v1", ""], 
#["", "e2-f.sub, e2-r.sub, cCcC", "v2, e3-f.sub, e3-r.sub, dDdD"]]
#adj=[[1, 1, 0],[0 ,0 ,0], [0 ,1 ,1]]
#vfiles, efiles, isoms, init_term_vs=unpack_gog(file_lines, adj)
#print(vfiles)
#print(efiles)
#print(isoms)
#print(init_term_vs)