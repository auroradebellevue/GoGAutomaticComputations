#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 11:49:45 2020

@author: aurora

This script contains functions which run kbmag on each rws ans sub rws in a 
particular directory. This script additionally contains functions to aid 
processing the input files.  

"""

##############################################################################
# imported files

import subprocess
import os

##############################################################################
# global variables
#kbmag_ftn_dir="/home/aurora/gap-4.11.1/pkg/kbmag-1.5.9/standalone/bin/"

#diry=sys.argv 
dir_path = os.path.dirname(os.path.realpath(__file__)) 

#############################################################################
# 
def find_files(dir_path, pattern, flag):
    """
    Parameters
    ----------
    dir_path: absolute directory path
    pattern: pattern to look for in the file
    flag: can be 'start' or 'end' to determine where to look for the pattern 
    in the file name

    Returns
    -------
    file_list: a list of files found

    """
    file_list=[]
    if flag=='start':
        for root, dirs, files in os.walk(dir_path): 
            for file in files:  
                if file.startswith(pattern):
                    file_list.append(root+str(file))
    elif flag=='end':
        for root, dirs, files in os.walk(dir_path):
            for file in files:
                if file.endswith(pattern):
                    file_list.append(root+str(file))
    else:
        file_list.append("Failed to find file!")        
    return file_list 

def basic_processing(file_object, flag):
    """
    Parameters
    ----------
    file_object : TYPE
        DESCRIPTION.

    Returns
    -------
    The lines of the file as a list of strings. Empty lines and commented lines
    have been deleted.

    """
    file_lines=open(file_object[0], 'r')
    #read the file in and split by row
    file_lines=file_lines.read().split('\n')
    #remove the lines of the header which start with #
    file_lines=list(filter(lambda x: (x.startswith('#') == False), file_lines))
    #remove empty lines
    file_lines=list(filter(lambda x: (x != ""), file_lines))
    if flag=='files':
        for i in range(len(file_lines)):
            file_lines[i]=file_lines[i].split('&')
    return file_lines

def get_graph_and_tree(graph_info):
    """
    This function takes the graph information as a list of strings. Each entry
    of the list represents a row of the adjacency matrix. 

    Parameters
    ----------
    graph_info : TYPE
        DESCRIPTION.

    Returns
    -------
    An adjacency matrix as a list of lists. 

    """
    for i in range(len(graph_info)-1):
        #split by a space which marks the different columns for the same row
        graph_info[i]=graph_info[i].split(' ')
        #print(graph_info[i])
        for j in range(len(graph_info[i])):
            graph_info[i][j]=int(graph_info[i][j])
    if graph_info[-1]=="No tree":
        tree=[]
    else:
        temp=graph_info[-1].replace(" ", "")
        tree=temp.split(',')
    return graph_info[:-1], tree
   
#############################################################################
# run autgroup on a files. note that the subgroup name must be the same as the 
# group name
def run_auts(gog, vfiles, efiles, init_term_vs, file_dir, kbmag_ftn_dir):
    
    for v in gog.vertices:
        subprocess.run([kbmag_ftn_dir+"autgroup", file_dir+v.group.file["vx"]])
        #TODO: run autgroup on the other orderings
    for e in gog.edges:
        #forward oriented subgroup
        print('\n', "Running autgroup then autcos on", file_dir + e.fg.supergp_file,'\n')
        subprocess.run([kbmag_ftn_dir+"autgroup", file_dir+e.fg.supergp_file])
        subprocess.run([kbmag_ftn_dir+"autcos", file_dir+e.fg.supergp_file])
        #reverse oriented subgroup
        print('\n', "Running autgroup then autcos on", file_dir + e.rg.supergp_file, '\n')
        subprocess.run([kbmag_ftn_dir+"autgroup", file_dir+e.rg.supergp_file])
        subprocess.run([kbmag_ftn_dir+"autcos", file_dir+e.rg.supergp_file])
    #run autgroup on the subgroup of the base vertex
    e0gp=gog.vertices[0].subgp
    print('\n', "Running autgroup ", file_dir, e0gp, '\n')
    subprocess.run([kbmag_ftn_dir+"autgroup", file_dir+e0gp])

##############################################################################
# testing ground
#gplist=find_files("/home/aurora/Documents/PythonGroupTheory/HNN-examples/")
#run_auts(gplist)

#s=subprocess.check_output(["find", "/home/aurora/Documents/PythonGroupTheory", "-type", "f","-name", "*.sub"])
#t=subprocess.check_output(["find", "/home/aurora/Documents/PythonGroupTheory", "-name", "vgp*"])
#print(s)
#print(t) 


