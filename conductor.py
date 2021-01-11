#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 16:11:59 2020

@author: Aurora Marks

This python script is the main file to "conduct" and control the processes of 
the files GraphofGroups.py, run-auto.py, HigginsNF.py, and FSA_manipulations.py.

There are three main stages of this procedure: 
    I. Manage the input and internal storage of the graph of groups
    II. Run the appropriate kbmag programs to obtain the automatic structures
    III. Form the Word Acceptor for Normal Forms of the Fundamental Group of
        the graph of groups. 
    IV. Form the Multiplier Automaton
    
    Notes on each stage:
        I. This python script has one input, the absolute path of the directory
        that contains all group rewriting files and subgroup rewriting files. 
        
"""
#import libraries
import sys

#import other scripts
#import GraphofGroups as GraphofGroups
import run_auto as runauto
#import HigginsNF as HigginsNF
#import FSA_manipulations as FSA 


def main(argv):
    """
    Parameters
    ----------
    argv : absolute path of the directory that contains the group files

    Returns
    -------
    Something evenutally

    """
    # Input and Internal Storage
    file_dir=argv[0]
    # Find file containing graph of groups structure
    graph_input=runauto.find_files(dir_path=file_dir, pattern='graph.txt', flag='start')
    graph_file=open(graph_input[0], 'r')
    #read the file in and split by row
    graph=graph_file.read().split('\n')
    CHANGE#remove the lines of the header which start with #
    CHANGE graph=list(filter(lambda x: (x.startswith('#') == False), graph))
    CHANGE#remove empty lines
    CHANGE graph=list(filter(lambda x: (x != ""), graph))
    CHANGE#change graph to the only the string containing the graph info
    CHANGE graph=graph[0:]
    print(graph)
    CHANGE for i in range(len(graph)):
        #split by '&' which marks the different information by column for the same row
        if '&' in graph[i]:
            graph[i]=graph[i].split('&')
        for j in range(len(graph[i])):
            #split by ',,' which marks the different information for the same row and column.
            graph[i]=graph[i].split(',,')
            #split by ',' which separates the subgroup file name and isomorphism description
            #graph[i][j][1]=graph[1].split(',')
    print(graph, type(graph), len(graph))
    
    # Run kbmag programs

    # Word Acceptor
    
if __name__=="__main__":
    main(sys.argv[1:]) #the first argument is the python script
