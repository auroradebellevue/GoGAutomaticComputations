#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 21:01:55 2020

@author: aurora
"""
import numpy as np

class Group():
    """This class holds the generating set and relator set for a group. It 
    should be designed to make everything easier for the user. In other words,
    a user can input a group presentation and not worry about inverses. 
    g_genset should be given as a list of lowercase letters, uppercase letters
    will be used for inverse letters.
    g_relset should be given as words equal to the identity"""
    def __init__(self, g_genset, g_relset, gen_set_order, order):
        self.g_genset=[]
        self.g_relset=[]
        self.gen_set_order=[]
        self.order='shortlex'
        self.initial_rules=[]
        
        for word in g_relset:
            self.initial_rules= np.append(self.initial_rules, [word, 1])
    
class Rules():
    """This class holds rules for an instance of a group. Each rule is stored
    as a row where the left entry is the left side of the rule and the right
    entry is the right side of the rule"""
    
G=Group(['a','b'],['abAB'],['a','b'], 'shortlex')
print(G)
print(G.initial_rules)    