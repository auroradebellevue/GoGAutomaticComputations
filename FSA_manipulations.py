#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 20:19:18 2020

@author: aurora

The functions in this script should be used to store FSA's in python and prep
the FSA's so that they're ready to be used by HigginsNF.py
"""
##############################################################################
# imported files
import subprocess
import os
##############################################################################
# global variables
kbmag_std_al_dir="/home/aurora/gap-4.11.0/pkg/kbmag-1.5.9/standalone/"
kbmag_ftn_dir="/home/aurora/gap-4.11.0/pkg/kbmag-1.5.9/standalone/bin/x86_64-pc-linux-gnu-default64-kv7/"
############################################################################
# class objects 
class fsa(object):
    """ 
    class to mimic an fsa!
    The alphabet should be an instance of the class alphabet. 
    The initial state should be a list. Unless the fsa is a mipda, this
    list will have one element. The attribute "accepting" should be a list
    of accepting states. The table should be an instance of the class table. 
    """
    def __init__(self, alphabet, initial, accepting, table):
        self.alphabet=alphabet 
        #self.states=states
        #self.flags TODO: there are many flags to watch for
        self.initial=initial
        self.accepting=accepting
        self.table=table
        #self.filename, TODO: creating a filename attribute seems like it'll come in handy soon :)
    
    def new_letter(self, l):
        """ This function adds a new letter and its inverse to the
        alphabet and no new transitions to the table"""
        self.alphabet.names.extend([f"{l}", f"{l.upper()}"])
        self.alphabet.size += 2
        for row in self.table.transitions:
            row.append(0)
            row.append(0)
    
    def permute_order(self, perm):
        """ 
        This function permutes the order of the alphabet and the columns
        of the transition table according the to permutation perm, represented
        as a list. The number in position 0 indicates the new index of the 
        first letter of the alphabet, etc.
        """
        #initialize temp_list with the correct length
        temp_list=[None]*len(perm)
        #reorder alphabet
        for i in range(len(perm)):
            temp_list[perm[i]]=self.alphabet.names[i]
        self.alphabet.names=temp_list
        #reorder transitions
        for i in range(len(self.table.transitions)):
            temp_list=[None]*len(perm)
            for j in range(len(perm)):
                temp_list[perm[j]]=self.table.transitions[i][j]
            self.table.transitions[i]=temp_list 
        #TODO: print out reordered fsa with an appropriate name
        #self.print_fsa("somefilename.wa")
            
    def print_fsa(self, filename):
        """
        This function prints out an fsa object to a file that kbmag can 
        interpret. The name filename should have the correct extension at the
        end when this method is called.

        """        
        short_tab="        " #8 spaces
        long_tab="                " #16 spaces
        #create file 
        file=open(filename, 'w')
        try:
            file.write("_RWS.wa := rec(\n")
            file.write(short_tab+"isFSA := true,\n")
            file.write(short_tab+"alphabet := rec(\n")
            file.write(long_tab+"type := \"identifiers\",\n")
            file.write(long_tab+"size := "+str(self.alphabet.size)+",\n")
            file.write(long_tab+"format := \"dense\",\n")
            self.alphabet.print_names(file, long_tab)
            file.write(long_tab+"),\n")
            file.write(short_tab+"states := rec(\n")
            file.write(long_tab+"type := \"simple\",\n")
            file.write(long_tab+"size := "+str(self.table.state_size)+"\n")
            file.write(long_tab+"),\n")
            file.write(short_tab+"flags := [\"DFA\",\"minimized\",\"BFS\",\"accessible\",\"trim\"],\n")
            file.write(short_tab+"initial := "+str(self.initial)+",\n")
            file.write(short_tab+"accepting := "+str(self.accepting)+",\n")
            file.write(short_tab+"table := rec(\n")
            file.write(short_tab+"format := \""+str(self.table.format)+"\",\n")
            file.write(short_tab+"numTransitions := "+str(self.table.num_transitions)+",\n")
            file.write(short_tab+"transitions := "+str(self.table.transitions))
            file.write(")\n")
            file.write(");")
        finally:
            file.close()
        
class alphabet(object):
    """class to hold all information about an alphabet!"""
    def __init__(self, size, names):
        self.size=size
        self.names=names #in order as a list
    
    def print_names(self, file, long_tab):
        """
        prints the names of the alphabet with the correct format to the file
        object, file.

        """
        file.write(long_tab+"names := [")
        for i in range(self.size):
            file.write(str(self.names[i]))
            if i!=self.size-1:
                file.write(",")
            else:
                file.write("]\n")
    def create_star_fsa(self):
        """
        creates an fsa that accepts all words over the alphabet.
        This fsa is saved to a file called "alphabet.wa.star"

        """
        t=[]
        accepting=[]
        temp_row=[]
        #the first row should describe the transitions from state 1 to a 
        #unique state for each letter of the alphabet, starting at state 2
        for i in range(len(self.names)):
            accepting.append(i+2)    
        t.append(accepting)
        for i in range(len(accepting)):
            temp_row.append(0)
        #We need rows of 0's for each state we created.
        for _ in temp_row:
            t.append(temp_row)
        full_table=table('dense deterministic', t)
        one_letter_fsa=fsa(self, [1], accepting, full_table)
        one_letter_fsa.print_fsa("alphabet.wa") #TODO we'll probably need to change this name later
        subprocess.run([kbmag_ftn_dir+"fsastar", "/home/aurora/Documents/PythonGroupTheory/alphabet.wa"],
                        stdout=subprocess.PIPE)
        
class table(object):
    """class to hold all information about the transition table!"""
    def __init__(self, table_format, transitions):
        self.format=table_format
        self.transitions=transitions
        self.state_size=len(self.transitions)
        self.num_transitions=self.count()
    
    def count(self):
        total=0
        for row in self.transitions:
            for entry in row:
                if entry!=0:
                    total+=1
        return total
##############################################################################
# Functions
def create_small_fsa(alphabet, letter):
    """
    creates an fsa that accepts a single letter in the alphabet given. The fsa
    will have one transition from state 1 to state 2, and state 2 will be the
    only accepting state. state 0 is the fail state.
    """
    if letter in alphabet.names:
        index=alphabet.names.index(letter)
        t=[]
        row0=[]
        row1=[]
        for i in range(alphabet.size):
            if i==index:
                row0.append(2)
                row1.append(0)
            else:
                row0.append(0)
                row1.append(0)
        t.append(row0)
        t.append(row1)
        full_table=table('dense deterministic', t)
        return fsa(alphabet, [1], [2], full_table)
    else:
        print("The alphabet does not contain this letter.")
        
def create_fsa_from_file(filename):
    """
    This function should read in a fsa file from kbmag and create an instance
    of an fsa object.
    """
    file=open(filename, "r")
    try:
        #read in all lines and store them in a list
        all_lines=file.readlines()
        
        #delete any lines that were comments
        index_list=[]
        for i in range(len(all_lines)):
            if all_lines[i].find("#")!=-1:
                index_list.append(i)
        index_list.sort(reverse=True)
        for i in index_list:
            del all_lines [i:i+1]
        
        #as implemented, the first 4 lines are not needed
        for _ in range(4):
            all_lines.pop(0)
        
        #retrieve size of alphabet, note the line ends with the 2 characters ,\n
        info_index=all_lines[0].find("=")
        alpha_size=int(all_lines[0][info_index+1:-2])
        all_lines.pop(0)
        all_lines.pop(0)
        
        #retrieve athe alphabet and store it as a list of strings
        info_index=all_lines[0].find("[")
        alpha_names=[]
        for i in range(alpha_size):
            alpha_names.append(all_lines[0][info_index+2*i+1:info_index+2*i+2])
        
        #create alphabet object
        alpha=alphabet(alpha_size, alpha_names)
        all_lines.pop(0)
        
        #as implemented, the next 3 lines are not needed
        for _ in range(3):
            all_lines.pop(0)
        
        #retrieve the number of states, note the line ends with the 1 character \n
        info_index=all_lines[0].find("=")
        fsa_state_size=int(all_lines[0][info_index+2:-1])
        all_lines.pop(0)
        
        #as implemented, the next 2 lines are not needed
        for _ in range(2):
            all_lines.pop(0)
        
        #retrieve initial state #TODO: this needs modification for mipda's
        info_index=all_lines[0].find("[")
        fsa_initial=[all_lines[0][info_index+1:-3]]
        fsa_initial[0]=int(fsa_initial[0])
        all_lines.pop(0)
        
        #retrieve accepting states
        info_index=all_lines[0].find("[")
        fsa_accepting=[]
        #if the accepting states are a list written in shorthand because the
        #accepting state labels are consecutive numbers
        if all_lines[0].find("..")!=-1:
            fsa_accepting.append(int(all_lines[0][info_index+1:info_index+2]))
            last_state=int(all_lines[0][info_index+4:info_index+5])
            i=int(all_lines[0][info_index+1:info_index+2])+1
            while i <=last_state:
                fsa_accepting.append(i)
                i+=1
        #if the accepting states are a list written with commas 
        elif all_lines[0].find("]")-info_index>2:
            temp_index=all_lines[0].find("]")
            temp_str=all_lines[0][info_index+1:temp_index]
            l=temp_str.split(",")
            for i in range(len(l)):
                l[i]=int(l[i])
            fsa_accepting=l
        #if there is one accept state
        else :   
            fsa_accepting.append(int(all_lines[0][info_index+1:info_index+2]))
        all_lines.pop(0)
        
        #as implemented, the next 3 lines are not needed (assume all transition
        #tables are dense deterministic)
        for _ in range(3):
            all_lines.pop(0)
        
        #retrieve transition table
        fsa_transitions=[]
        for i in range(len(all_lines)):
            all_lines[i]=all_lines[i].strip()
        #join all lines of the code just so that the whole table is in one entry
        str_table=' '.join(all_lines[:])
        #remove spaces
        str_table=str_table.replace(" ", "")
        #remove the words, :=, extra ));, and [ ]
        str_table=str_table[13:-3]
        str_table=str_table.replace("[", "")
        str_table=str_table.replace("]", "")
        #split str_table by commas to give list as strings
        str_table=str_table.split(",")
        #join together entries according to the size of the alphabet and 
        #convert entries from strings to integers (for printing later)
        for i in range(fsa_state_size):
             temp_row=str_table[0+i*alpha_size : alpha_size*(i+1)]
             for j in range(len(temp_row)):
                 temp_row[j]=int(temp_row[j])
             fsa_transitions.append(temp_row)
        
        #create table object
        t=table("dense deterministic", fsa_transitions)
        
        #create fsa object
        return fsa(alpha, fsa_initial, fsa_accepting, t)
    finally:
        file.close()
        
#TODO: create function that compares alphabets and determines the permutation
        #needed to get the ordering to match the 2nd one to the 1st one.
"""
##############################################################################
# Driver code for testing
        
a=alphabet(4, ['a', 'A', 'b', 'B'])
t=table('dense deterministic', [[2,3,4,5],
                         [2,0,4,5],
                         [0,3,4,5],
                         [0,0,4,0],
                         [0,0,0,5] 
                        ]) 
f=fsa(alphabet=a, initial=[1],accepting=[1, 2, 3, 4, 5], table=t)

f.new_letter('t')
#f.permute_order([1,0,2,3,4,5])
#print(f.alphabet.names)
#print(f.table.transitions)

newfsa=create_small_fsa(a, 'A')
#print(newfsa.table.transitions)
newfsa.print_fsa("newfsa.wa")

#a.create_star_fsa()
#fsa_alphabet=create_fsa_from_file("alphabet.wa")
#print(fsa_alphabet)
##############################################################################
# testing area
#trial=subprocess.run(["/home/aurora/gap-4.11.0/pkg/kbmag-1.5.9/standalone/bin/x86_64-pc-linux-gnu-default64-kv7/autgroup", "/home/aurora/gap-4.11.0/pkg/kbmag-1.5.9/standalone/ag_data/236"],
#                     stdout=subprocess.PIPE)
#newtrial=subprocess.run([kbmag_ftn_dir+"fsaconcat", "/home/aurora/Documents/PythonGroupTheory/newfsa.wa", "/home/aurora/Documents/PythonGroupTheory/newfsa.wa", "/home/aurora/Documents/PythonGroupTheory/concat-test.wa"],
#                        stdout=subprocess.PIPE)
#print(newtrial)

test=create_fsa_from_file("Z2_bFirst.cos.wa")
test.new_letter('t')
test.print_fsa("Z2_bFirst_full_alpha.cos.wa")
"""