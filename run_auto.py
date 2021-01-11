#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 11:49:45 2020

@author: aurora
THESE FUNCTIONS NEED TO BE REPURPOSED FOR NEW ORGANIZATION NOW THAT GRAPH.TXT
WILL BE USED TO STORE THE VERTEX GROUPS, EDGE GROUPS, AND GRAPH INFO (
WITH NAMES IN THE GRAPH.TXT FILE)
This file runs kbmag on each rws ans sub rws in a folder. 

"""

##############################################################################
# imported files

import subprocess
import os
import sys

##############################################################################
# global variables
kbmag_std_al_dir="/home/aurora/gap-4.11.0/pkg/kbmag-1.5.9/standalone/"
kbmag_ftn_dir="/home/aurora/gap-4.11.0/pkg/kbmag-1.5.9/standalone/bin/x86_64-pc-linux-gnu-default64-kv7/"

#diry=sys.argv 
dir_path = os.path.dirname(os.path.realpath(__file__)) 

#############################################################################
# 
def OLD_find_files(diry):
    """
    assign the directory to diry and find vertex group files and edge group files
    then assign these file names to different lists

    Parameters
    ----------
    diry : TYPE
        DESCRIPTION.

    Returns
    -------
    A list with 2 entries. The first entry is the list of vertex group files. 
    The second entry is the list of edge group files. 

    """
    #s=subprocess.check_output(["find", diry, "-type", "f","-name", "*.sub"])
    #print(s)
    voutput=subprocess.check_output(["find", diry, "-name", "vgp\\*"])
    eoutput=subprocess.check_output(["find", diry, "-type", "f","-name", "*.sub"])
    return [voutput, eoutput]



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

#############################################################################
# run autgroup on a files. note that the subgroup name must be the same as the 
# group name
def run_auts(gplist):
    for file in gplist[0]:
        subprocess.run([kbmag_ftn_dir+"autgroup", file])
    for file in gplist[0]:
        subprocess.run([kbmag_ftn_dir+"autcos", file])

##############################################################################
# testing ground
#gplist=find_files("/home/aurora/Documents/PythonGroupTheory/HNN-examples/")
#run_auts(gplist)

#s=subprocess.check_output(["find", "/home/aurora/Documents/PythonGroupTheory", "-type", "f","-name", "*.sub"])
#t=subprocess.check_output(["find", "/home/aurora/Documents/PythonGroupTheory", "-name", "vgp*"])
#print(s)
#print(t)