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

import subprocess
import os
import FSA_manipulations as fsa_ftns

##############################################################################
# global variables
kbmag_std_al_dir="/home/aurora/gap-4.11.0/pkg/kbmag-1.5.9/standalone/"
kbmag_ftn_dir="/home/aurora/gap-4.11.0/pkg/kbmag-1.5.9/standalone/bin/x86_64-pc-linux-gnu-default64-kv7/"
#############################################################################
#TODO:This process needs to be automated 
#TODO:Need a better naming scheme 
#one vertex, one edge
#\[\hat{N}= L(WaR)\cdot \left( (s_e L(WaS1)) \cup (s_{\overline{e}} L(WaS2) )\right)^{*} \setminus \left[ \left(\hat{X}^{*} s_{e}s_{\overline{e}} \hat{X}^{*}\right) \cup \left(\hat{X}^{*} s_{\overline{e}}s_e \hat{X}^{*} \right)\right] \]

#read in vertex group word acceptor and edge group coset word acceptors
wa_vertex=fsa_ftns.create_fsa_from_file("Z2_bFirst.wa")
wa_edge=fsa_ftns.create_fsa_from_file("Z2_bFirst.cos.wa")
wa_edge_rev=fsa_ftns.create_fsa_from_file("Z2_aFirst.cos.wa")

#reorder reverse edge to match alphabets
wa_edge_rev.permute_order([2,3,0,1])

#add stable letter with inverse to the end of the alphabet of each of these fsa's
wa_vertex.new_letter('t')
wa_edge.new_letter('t')
wa_edge_rev.new_letter('t')

#print to file the fsa's with their full alphabet
wa_vertex.print_fsa("Z2_bFirst_fullalpha.wa")
wa_edge.print_fsa("Z2_bFirst_fullalpha.cos.wa")
wa_edge_rev.print_fsa("Z2_aFirst_fullalpha.cos.wa")

#create small fsa's needed for normal forms and print to file
stab_wa=fsa_ftns.create_small_fsa(wa_vertex.alphabet, 't')
stab_wa.print_fsa("stab.wa")
stab_inv_wa=fsa_ftns.create_small_fsa(wa_vertex.alphabet, 'T')
stab_inv_wa.print_fsa("stab_inv.wa")
subprocess.run([kbmag_ftn_dir+"fsaconcat", 
                "/home/aurora/Documents/PythonGroupTheory/stab.wa",
                "/home/aurora/Documents/PythonGroupTheory/stab_inv.wa",
                "/home/aurora/Documents/PythonGroupTheory/tT.wa"])
subprocess.run([kbmag_ftn_dir+"fsaconcat",
                "/home/aurora/Documents/PythonGroupTheory/stab_inv.wa",
                "/home/aurora/Documents/PythonGroupTheory/stab.wa",
                "/home/aurora/Documents/PythonGroupTheory/Tt.wa"])
wa_vertex.alphabet.create_star_fsa()

#use the functions in kbmag to create the normal forms
##############################################################################
#product of stable letter with coset reps
print("t (b* or B*)")
subprocess.run([kbmag_ftn_dir+"fsaconcat",
                "/home/aurora/Documents/PythonGroupTheory/stab.wa",
                "/home/aurora/Documents/PythonGroupTheory/Z2_aFirst_fullalpha.cos.wa",
                "/home/aurora/Documents/PythonGroupTheory/stab_t_with_cosets.wa"])
print("T (a* or A*)")
subprocess.run([kbmag_ftn_dir+"fsaconcat",
                "/home/aurora/Documents/PythonGroupTheory/stab_inv.wa",
                "/home/aurora/Documents/PythonGroupTheory/Z2_bFirst_fullalpha.cos.wa",
                "/home/aurora/Documents/PythonGroupTheory/stab_T_with_cosets.wa"])
print("t(b* or B*) OR T(a* or A*)")
subprocess.run([kbmag_ftn_dir+"fsaor",
                "/home/aurora/Documents/PythonGroupTheory/stab_t_with_cosets.wa",
                "/home/aurora/Documents/PythonGroupTheory/stab_T_with_cosets.wa",
                "/home/aurora/Documents/PythonGroupTheory/stab_letters_with_cosets.wa"])
print("(t(b* or B*) OR T(a* or A*))*")
subprocess.run([kbmag_ftn_dir+"fsastar",
                "/home/aurora/Documents/PythonGroupTheory/stab_letters_with_cosets.wa"])
###############################################################################
#Building bad words
print("X* tT")
subprocess.run([kbmag_ftn_dir+"fsaconcat",
                "/home/aurora/Documents/PythonGroupTheory/alphabet.wa.star",
                "/home/aurora/Documents/PythonGroupTheory/tT.wa",
                "/home/aurora/Documents/PythonGroupTheory/anyword_tT.wa"])
print("X* tT X*")
subprocess.run([kbmag_ftn_dir+"fsaconcat",
                "/home/aurora/Documents/PythonGroupTheory/anyword_tT.wa",
                "/home/aurora/Documents/PythonGroupTheory/alphabet.wa.star",
                "/home/aurora/Documents/PythonGroupTheory/anyword_tT_anyword.wa"])
print("X* Tt")
subprocess.run([kbmag_ftn_dir+"fsaconcat",
                "/home/aurora/Documents/PythonGroupTheory/alphabet.wa.star",
                "/home/aurora/Documents/PythonGroupTheory/Tt.wa",
                "/home/aurora/Documents/PythonGroupTheory/anyword_Tt.wa"])
print("X* Tt X*")
subprocess.run([kbmag_ftn_dir+"fsaconcat",
                "/home/aurora/Documents/PythonGroupTheory/anyword_Tt.wa",
                "/home/aurora/Documents/PythonGroupTheory/alphabet.wa.star",
                "/home/aurora/Documents/PythonGroupTheory/anyword_Tt_anyword.wa"])
print("X* tT X* OR X* Tt X")
subprocess.run([kbmag_ftn_dir+"fsaor",
                "/home/aurora/Documents/PythonGroupTheory/anyword_tT_anyword.wa",
                "/home/aurora/Documents/PythonGroupTheory/anyword_Tt_anyword.wa",
                "/home/aurora/Documents/PythonGroupTheory/anyword_tT_or_Tt_anyword.wa"])
###############################################################################
# Final Steps
subprocess.run([kbmag_ftn_dir+"fsaconcat", 
                "/home/aurora/Documents/PythonGroupTheory/Z2_bFirst_fullalpha.wa",
                "/home/aurora/Documents/PythonGroupTheory/stab_letters_with_cosets.wa.star",
                "/home/aurora/Documents/PythonGroupTheory/vertex_words_with_stab_letters_and_cosets.wa"])

subprocess.run([kbmag_ftn_dir+"fsaandnot",
                "/home/aurora/Documents/PythonGroupTheory/vertex_words_with_stab_letters_and_cosets.wa",
                "/home/aurora/Documents/PythonGroupTheory/anyword_tT_or_Tt_anyword.wa",
                "/home/aurora/Documents/PythonGroupTheory/HigginsNF.wa"])


