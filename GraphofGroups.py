#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 16:31:45 2020

@author: aurora
"""
class GoG(object):
    """
    This is a graph of groups.
    """
    def __init__(self, v, e, i):
        self.vertices=[v] #TODO: use input list for vertex set
        self.edges=[e] #TODO: use list-of-lists to store edge info 
        self.incidence=i #TODO: use an incidence matrix to keep track of graph shape(?)

class vertex(object):
    """
    This is a vertex of a graph of groups.
    """
    def __init__(self, label, G):
        self.label=label
        self.group=G

class edge(object):
    """
    This is a edge and the reverse edge of a graph of groups.
    """
    def __init__(self, label, iv, tv, fG, rG):
        self.label=label
        self.initial_vertex=iv
        self.terminal_vertex=tv
        self.forward_group=fG
        self.reverse_group=rG

class group(object):
    """
    This is a group.
    """
    def __init__(self, gens, invgens, rels):
        self.gens=gens
        self.invgens=invgens
        self.rels=rels
        self.ordering="shortlex"
        
#############################################################################
# Testing area/driver code
Gv=group(['a', 'A', 'b', 'B'], ['A', 'a', 'B', 'b'], ['abAB'])
v=vertex(1, Gv)
Gfe=group(['a', 'A'], ['A', 'a'], [])
Gre=group(['b', 'B'], ['B', 'b'], [])
e=edge(1, v, v, Gfe, Gre)
i=[]
gog=GoG(v, e, i)
