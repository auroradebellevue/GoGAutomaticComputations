#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 20:19:18 2020

@author: Aurora Marks and Ash DeClerk

The functions in this script should be used to store FSA's in python and prep
the FSA's so that they're ready to be used by HigginsNF.py
"""
##############################################################################
# imported files
import subprocess
import re
##############################################################################
# global variables
#kbmag_ftn_dir="/home/aurora/gap-4.11.1/pkg/kbmag-1.5.9/standalone/bin/"
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
    def __init__(self, alpha, sts, tbl):
        self.alphabet=alpha
        self.states=sts
        self.flags=[]
        self.table=tbl
        self.file_lines=""
        self.filename="" 

    def new_letter(self, l):
        """ This function adds a new letter to the alphabet and no new 
        transitions to the table. The inverse must be added manually."""
        self.alphabet.names.extend(l) 
        self.alphabet.size += 1
        for row in self.table.transitions:
            row.append(0)
    
    def new_stab_letter(self, l):
        """ This function adds a new stable letter and its inverse to the
        alphabet and no new transitions to the table. The stable letter
        should look like t1, t2, t(number). """
        label=l[1:]
        self.alphabet.names.extend(["t"+label, "T"+label]) 
        self.alphabet.size += 2
        for row in self.table.transitions:
            row.append(0)
            row.append(0)
            
            
    def print_fsa(self, filename):
        """
        This function prints out an fsa object to a file that kbmag can 
        interpret. The name filename should have the correct extension at the
        end when this method is called.

        """     
        sTab="        " #8 spaces
        #create file 
        file=open(filename, 'w')
        try:
            if "wa" in self.flags:
                file.write("_RWS.wa := rec(\n")
            elif "gm" in self.flags:
                file.write("_RWS.gm := rec(\n")
            else: 
                file.write("_RWS.mult := rec(\n")
            file.write(sTab+"isFSA := true,\n")
            if self.alphabet.type=="product":
                file.write(sTab+"alphabet := rec(\n")
                file.write(sTab+sTab+"type := \""+ str(self.alphabet.type)+"\", \n")
                file.write(sTab+sTab+"size := "+str(self.alphabet.size)+",\n")
                file.write(sTab+sTab+"arity := 2,\n")
                file.write(sTab+sTab+"padding := _,\n")
                file.write(sTab+sTab+"base := rec( \n")
                file.write(sTab+sTab+"type := \"identifiers\",\n")
                file.write(sTab+sTab+"size := "+ str(len(self.alphabet.base))+", \n")
                file.write(sTab+sTab+"format := \"dense\",\n")
                file.write(sTab+sTab+"names := [")
                self.alphabet.print_names(file, sTab+sTab, 'product')
                file.write(sTab+sTab+") \n" )
            else:#all other cases    
                file.write(sTab+"alphabet := rec(\n")
                file.write(sTab+sTab+"type := \""+ str(self.alphabet.type)+"\", \n")
                file.write(sTab+sTab+"size := "+str(self.alphabet.size)+",\n")
                file.write(sTab+sTab+"format := \"dense\",\n")
                file.write(sTab+sTab+"names := [")
                self.alphabet.print_names(file, sTab+sTab, 'dense')
            file.write(sTab+sTab+"),\n")
            if self.states.format !='labeled': #no labeled states
                file.write(sTab+"states := rec(\n")
                file.write(sTab+sTab+"type := \"simple\",\n")
                file.write(sTab+sTab+"size := "+str(self.states.size)+"\n")
                file.write(sTab+sTab+"),\n")
                file.write(sTab+"flags := [\"DFA\"],\n")
                file.write(sTab+"initial := ["+str(self.states.initial[0])+"],\n")
                file.write(sTab+"accepting := "+str(self.states.accepting)+",\n")
            else: #states are labeled
                file.write(sTab+"states := rec(\n")
                file.write(sTab+ "type :=\"labeled\",\n")
                file.write(sTab+ "size := "+ str(self.states.size)+",\n")
                if len(self.states.label_lines)!=0:
                    for line in self.states.label_lines:
                        line=line.replace("\n", "")
                        file.write(line+ "\n")
                else:
                    file.write(sTab + "labels := rec( \n")
                    file.write(sTab +sTab + "type := \"list of words\",\n")
                    file.write(sTab +sTab + "size := " + str(len(self.states.label_states_pairs))+",\n")
                    #print alphabet of base correctly
                    file.write(sTab +sTab + "alphabet := [")
                    size=len(self.alphabet.base)
                    for i in range(size):
                        file.write(str(self.alphabet.base[i]))
                        if i!=size-1:
                            file.write(",")
                        else:
                            file.write("],\n")
                    file.write(sTab + sTab + "format := \"sparse\",\n")
                    
                    #print integer labels of word difference labels
                    file.write(sTab+sTab+ "names := [ \n")
                    st_label_nums=len(self.states.label_numbers_pairs)
                    for i in range(st_label_nums):
                        file.write(sTab+sTab + '['+ str(self.states.label_numbers_pairs[i][0])
                                   +',['+str(self.states.label_numbers_pairs[i][1])+']]')
                        if i!= st_label_nums -1:
                            file.write(',\n')
                        else:
                            file.write('\n')
                    file.write(sTab+sTab+']\n')
                    file.write(sTab+sTab+"),\n")
                    
                    file.write(sTab+sTab+"format := \"sparse\",\n")
                    file.write(sTab+sTab+"setToLabels := [ \n")
                    st_label_pairs = len(self.states.label_states_pairs)
                    for i in range(st_label_pairs):
                        file.write(sTab+sTab + '['+ str(self.states.label_states_pairs[i][0])
                                   +','+str(self.states.label_states_pairs[i][1])+']')
                        if i!= st_label_pairs -1:
                            file.write(',\n')
                        else:
                            file.write('\n')
                    file.write(sTab+sTab+']\n')
                    file.write(sTab+sTab+"),\n")
                               
                file.write(sTab+"flags := [\"DFA\"],\n")
                file.write(sTab+"initial := ["+str(self.states.initial[0])+"],\n")
                file.write(sTab+"accepting := "+str(self.states.accepting)+",\n")
            if self.table.format=="dense deterministic":
                file.write(sTab+"table := rec(\n")
                file.write(sTab+"format := \""+str(self.table.format)+"\",\n")
                file.write(sTab+"numTransitions := "+str(self.table.num_transitions)+",\n")
                file.write(sTab+"transitions := \n")
                for i in range(len(self.table.transitions)):
                    if i==0:
                        file.write("[")
                    file.write(sTab+ str(self.table.transitions[i]))
                    if i!=len(self.table.transitions)-1:
                        file.write(",\n")
                file.write("] \n")
                file.write(")\n") #close table record
                file.write(");") #close fsa record
            else: #table is sparse
                file.write(sTab+"table := rec(\n")
                file.write(sTab+"format := \""+str(self.table.format)+"\",\n")
                file.write(sTab+"numTransitions := "+str(self.table.num_transitions)+",\n")
                file.write(sTab+"transitions := "+str(self.table.transitions) )
                file.write("\n ) \n")
                file.write(");")
        finally:
            file.close()

class ndfsa(object):
    """ 
    class to mimic an non deterministic fsa!
    The alphabet should be an instance of the class alphabet. 
    The initial state should be a list with one element. 
    The attribute "accepting" should be a list
    of accepting states. 
    The table should be an instance of the class ndtable. 
    """
    def __init__(self, alpha, sts, tbl):
        self.alphabet=alpha
        self.states=sts
        self.flags=[]
        self.table=tbl
        self.file_lines=""
        self.filename="" 
        self.epsilon=self.e_free_test()
    def e_free_test(self):
        """Tests if epsilon represented by _ is in the alphabet"""
        if '_' in self.alphabet.names:
            return False
        else:
            return True
class alphabet(object):
    """class to hold all information about an alphabet!"""
    def __init__(self, size, names, alpha_type):
        self.size=size
        self.names=names #in order as a list
        self.type=alpha_type #identifiers or product
        self.base=[] #used for product alphabets only
        
    def print_names(self, file, long_tab, style):
        """
        prints the names of the alphabet with the correct format to the file
        object, file.
        This function prints the full alphabet if it is dense and only prints
        the base list if it is a product alphabet.

        """
        if style == 'dense':
            for i in range(self.size):
                file.write(str(self.names[i]))
                if i!=self.size-1:
                    file.write(",")
                else:
                    file.write("]\n")
        elif style == 'product':
            size=len(self.base)
            for i in range(size):
                file.write(str(self.base[i]))
                if i!=size-1:
                    file.write(",")
                else:
                    file.write("]\n")
        
    def create_star_fsa(self, name, file_dir, kbmag_ftn_dir):
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
        sts=states(len(accepting)+1 , [1], accepting)
        one_letter_fsa=fsa(self, sts, full_table)
        one_letter_fsa.print_fsa(file_dir + name) 
        subprocess.run([kbmag_ftn_dir+"fsastar", file_dir + name],
                        stdout=subprocess.PIPE)
    def free_reduction(self, word):
        """
        This function outputs the freely reduced form of the input, word. 
        Free reduction turns pairs of letters "xX" into the empty string or
        "IdWord". 
        This function will only work for dense alphabets.

        Parameters
        ----------
        word : string that represents a word which needs to be freely reduced

        Returns
        -------
        reduced_word : the freely reduced word

        """
        if self.type == "product":
            return "ERROR; this function cannot perform free reduction over a product alphabet!"
        #different kinds of identity representations
        if word == "IdWord":
            return str()
        if len(word) == 0:
            return str()
        #else we may look for "lL" pairs
        pairs=[]        
        for l in self.names:
            pairs.append(l+l.swapcase())
        #print(pairs)
        #print("we word to simplify the word", word)
        IsFreelyReduced = [x in word for x in pairs]            
        #print(IsFreelyReduced)
        new_word = word
        while True in IsFreelyReduced:
            index1=IsFreelyReduced.index(True) #gives the first index 
            new_word = new_word.replace(pairs[index1], str()) #replace all occurances of lL with the empty string
            #print("Replaced: ", pairs[index1], " and obtained: ", new_word)
            if len(new_word) == 0:
                return str()
            IsFreelyReduced = [x in new_word for x in pairs]            
            
            #print("IsFreelyReduced", IsFreelyReduced)
            
        return new_word

class states(object):
    """
        class to hold the states, their labels, and other unique aspects
    """
    def __init__(self, size, init, accepting):
        self.size=size
        self.initial=init
        self.accepting=accepting
        self.format=""
        self.label_numbers_pairs = []
        self.label_states_pairs = []
        self.label_lines = ""
        
class table(object):
    """class to hold all information about the transition table!"""
    def __init__(self, table_format, transitions):
        self.format=table_format #dense deterministic or sparse
        self.transitions=transitions
        self.num_transitions=self.count()
    
    def count(self):
        total=0
        for row in self.transitions:
            for entry in row:
                if entry!=0:
                    total+=1
        return total
    
class ndtable(object):
    """class to hold all information about a non deterministic transition table!
    The data type for the table will need to use sets instead of the current format.
    It should be a 2D list of sets. The row represents the starting state, and 
    the column represents the transition letter and a set will be in that 
    position of the 2D list.
    """
    def __init__(self, table_format, transitions):
        self.format=table_format
        self.transitions=transitions
    
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
        sts=states(2, [1], [2])
        temp_fsa=fsa(alphabet, sts, full_table)
        return temp_fsa
    else:
        print("The alphabet does not contain this letter.")
def create_sub_fsa(M, alpha):
    """
    creates an fsa that accepts all words of length up to M over the alphabet
    alpha.
    """
    t=[]
    count=2
    for i in range(M):
        row=[]
        for j in range(len(alpha)):
            row.append(count)
        t.append(row)
        count+=1
    a=alphabet(len(alpha), alpha, "identifiers")
    tbl=table('dense deterministic', t)
    acc=list(range(1,M+2))
    sts=states(M, [1], acc)
    temp_fsa=fsa(a, sts,tbl)
    temp_fsa.flags=["dense deterministic"]
    return temp_fsa
    
def create_fsa_from_file(filename):
    """ 
    This function should read in a fsa file from kbmag and create an instance
    of an fsa object. It is a rewritten version of the old create_fsa function.
    This new function should work with more kinds of fsa files, especially
    files that have product alphabets or sparse tables. 

    """
    file=open(filename, "r")
    try:
        #read in all lines and store them in a list
        file_lines=file.readlines()

        #creates the alphabet!
        alpha=create_alphabet(file_lines)
        
        #creates the states and labels if applicable
        sts=create_states(file_lines, len(alpha.base))
        
        #creates the table and stores the type
        tbl=create_table(file_lines, alpha.size, sts.size)
        
        #create fsa object
        temp_fsa=fsa(alpha, sts, tbl)
        temp_fsa.file_lines=file_lines
        temp_fsa.filename=filename
        
        return temp_fsa
    finally:
        file.close()
def create_alphabet(lines):
    """
    Creates the alphabet of an fsa from the list of all lines from the file.
    This function should be used with the function "new_create_fsa_from_file"

    """
    
    #find sublist of lines that contain information about the alphabet
    for l in lines:
        if "alphabet" in l:
            index_1=lines.index(l,0)
            break
    for l in lines:
        if "states" in l:
            index_2=lines.index(l,0)
    alpha_lines=lines[index_1:index_2]
    #find the index of the lines for the type, size, and identifiers
    for l in alpha_lines:
        if "type" in l:
            type_index=alpha_lines.index(l) #first type distinguishes dense or product
            break
    for l in alpha_lines:
        if "size" in l:
            size_index=alpha_lines.index(l)
            break
    """
    for l in alpha_lines: 
        if "format" in l: 
            structure_index=alpha_lines.index(l)
    """
    for l in alpha_lines:
        if "names" in l:
            names_index=alpha_lines.index(l)
        if "simple" in l:
            names_index = None
    #find string description of the type
    info_index1=alpha_lines[type_index].find("\"")
    info_index2=alpha_lines[type_index].find("\"", info_index1+1)
    temp_string=alpha_lines[type_index][info_index1+1:info_index2]
    alpha_type=temp_string
    
    #find string description of the size
    info_index=alpha_lines[size_index].find("=")
    alpha_size=int(alpha_lines[size_index][info_index+1:].replace(",", ""))
    
    """
    #find string descriptions of the format/structure
    info_index1=alpha_lines[structure_index].find("\"")
    info_index2=alpha_lines[structure_index].find("\"", info_index1+1)
    temp_string=alpha_lines[structure_index][info_index1+1:info_index2]
    alpha_structure=temp_string
    """
    #find string list of the names
    if alpha_type != "simple":
        info_index1=alpha_lines[names_index].find("[")
        info_index2=alpha_lines[names_index].find("]")
        temp_string=alpha_lines[names_index][info_index1+1:info_index2]
        alpha_names=temp_string.split(",")
    
    #further processing for multplier automata
    if alpha_type=="product":
        temp_alpha=[]
        for l in alpha_names:
            for j in alpha_names:
                temp_alpha.append((l,j))
            temp_alpha.append((l, "_"))
        for l in alpha_names: 
            temp_alpha.append(("_", l))
        alpha_full_names=temp_alpha
    elif alpha_type =="simple":
        alpha_full_names = range(alpha_size)
    else:
        alpha_full_names=alpha_names
    a=alphabet(alpha_size, alpha_full_names, alpha_type)
    if alpha_type=="product":
        a.base=alpha_names
    return a
    
def create_states(lines, alpha_size):
    """
    Creates the table of an fsa from the list of all lines from the file.
    This function should be used with the function "new_create_fsa_from_file"

    """ 
    IsLabeled=False    
    #find sublist of lines that contain information about the states
    for l in lines:
        if "states" in l:
            index_1=lines.index(l,0)
            break
    for l in lines:
        if "table" in l:
            index_2=lines.index(l,0)
    state_lines=lines[index_1:index_2]
    #find and store the type of the state set
    for l in state_lines:
        if "type" in l:
            type_index=state_lines.index(l)
            break
    info_index1=state_lines[type_index].find("\"")
    info_index2=state_lines[type_index].find("\"", info_index1+1)
    temp_string=state_lines[type_index][info_index1+1:info_index2]
    state_structure=temp_string
    #find and store the size of the state set
    for l in state_lines:
        if "size" in l: #the first line with size is the state size
            state_size_index=state_lines.index(l)
            break
    info_index=state_lines[state_size_index].find("=")
    state_size=int(state_lines[state_size_index][info_index+1:].replace(",", ""))
    #find and store the initial state
    for l in state_lines:
        if "initial" in l: 
            init_state_index=state_lines.index(l)
    info_index1=state_lines[init_state_index].find("[")
    info_index2=state_lines[init_state_index].find("]")
    init_state=state_lines[init_state_index][info_index1+1:info_index2]
    init_state=[int(init_state)]
    #find and store the accepting states
    for l in state_lines:
        if "accepting" in l:
            accept_state_index=state_lines.index(l)
            accept_state_lines = state_lines[accept_state_index:]
            accept_state_lines = "".join(accept_state_lines)
    info_index1=accept_state_lines.find("[")
    info_index2=accept_state_lines.find("]")
    fsa_accepting=[]
    #if the accepting states are a list written in shorthand because the
    #accepting state labels are consecutive numbers
    if state_lines[accept_state_index].find("..")!=-1:
        temp_index=state_lines[accept_state_index].find("..")
        first_state=int(state_lines[accept_state_index][info_index1+1:temp_index])
        last_state=int(state_lines[accept_state_index][temp_index+2:info_index2])
        fsa_accepting=list(range(first_state, last_state+1))
    #if the accepting states are a list written with commas 
    elif state_lines[accept_state_index][info_index1+1:info_index2].find(",")!= -1:
        temp_str=accept_state_lines[info_index1+1:info_index2]
        print(temp_str)
        lst=temp_str.split(",")
        for i in range(len(lst)):
            lst[i]=int(lst[i])
        fsa_accepting=lst
    #if there is one accept state
    else :   
        fsa_accepting.append(int(state_lines[accept_state_index][info_index1+1:info_index2])) 
    #find the index of the lines for the label names and state labels
    for l in state_lines:
        if "labels" in l:
            labels_index=state_lines.index(l)
            IsLabeled=True
            break
    if IsLabeled==True:
        label_lines=state_lines[labels_index:-1]
        for l in label_lines:
            if "size" in l: #find index of size of labels
                info_index=l.find("=")
                label_size=int(l[info_index+1:].replace(",", ""))
                print("label size", label_size)
            if "names" in l: #find indices of list of label names
                label_names_index1=label_lines.index(l) + 1
                label_names_index2=label_names_index1+ alpha_size + 1
            if "setToLabels" in l: #find indices of label and state relationship
                label_scheme_index1=label_lines.index(l)
                label_scheme_index2=label_scheme_index1+len(fsa_accepting) + 1
        #create list of label names label i is in index i-1
        label_names=[]
        for l in label_lines[label_names_index1:label_names_index2]:
            first_brac=l.find("[")
            sec_brac=l.find("[", first_brac+1)
            end_brac=l.find("]")
            temp_list=[]
            temp_list.append(l[first_brac+1: sec_brac-1].replace(",", "")) 
            temp_list.append(l[sec_brac+1: end_brac].split(","))
            label_names.append(temp_list)
        #create list of labels for states
        state_labels=[]
        for l in label_lines[label_scheme_index1+1:label_scheme_index2]:
            beg_brac=l.find("[")
            end_brac=l.find("]")
            state_labels.append(l[beg_brac+1:end_brac].split(","))
            #size, init, accepting, label_names, state_labels)
        s=states(state_size, init_state, fsa_accepting)
        s.format=state_structure
        s.label_numbers_pairs=label_names
        s.label_states_pairs=state_labels
        s.label_lines=label_lines
    else:
        s=states(state_size, init_state, fsa_accepting)
    return s
    
def create_table(lines, alpha_size, state_size):
    """
    Creates the table of an fsa from the list of all lines from the file.
    This function should be used with the function "new_create_fsa_from_file"

    """ 
    #find sublist of lines that contain information about the states
    for l in lines:
        if "table" in l:
            index_1=lines.index(l,0)
            break
    table_lines=lines[index_1:]
    #find table format and store it as a string
    if table_lines[1].find("dense deterministic")!=-1:
        table_format="dense deterministic"
    elif table_lines[1].find("sparse")!=-1: 
        table_format="sparse"    
    else: 
        print("In the function create_table, table type is not supported")
    #find the transition table
    for l in table_lines: 
        if "transitions" in l:
            transition_index=table_lines.index(l)
    transition_lines=table_lines[transition_index: -1]
    if table_format=="dense deterministic":
        #retrieve transition table
        fsa_transitions=[]
        for i in range(len(transition_lines)):
            transition_lines[i]=transition_lines[i].strip() #delete leading and trailing space
        str_table=' '.join(transition_lines[:])#join all lines
        #remove inner spaces
        str_table=str_table.replace(" ", "")
        #remove the words, :=, extra ));, and [ ]
        bad_chars=[":", "=", "[", "]", ")", ";", "transitions"]
        for s in bad_chars:
            str_table=str_table.replace(s, "")
        #split str_table by commas to give list as strings
        str_table=str_table.split(",")
        #join together entries according to the size of the alphabet and 
        #convert entries from strings to integers (for printing later)
        for i in range(state_size):
            temp_row=str_table[0+i*alpha_size : alpha_size*(i+1)]
            for j in range(len(temp_row)):
                temp_row[j]=int(temp_row[j])
            fsa_transitions.append(temp_row)
        #create table object
        t=table("dense deterministic", fsa_transitions)

    elif table_format=="sparse":
        #retrieve transition table
        fsa_transitions=[]
        for i in range(len(transition_lines)):
            #delete leading and trailing space
            transition_lines[i]=transition_lines[i].strip() 
        #remove the words, :=, extra ));
        bad_chars=[":", "=", ")", ";", "transitions"]
        for i in range(len(transition_lines)):
            for s in bad_chars:
                transition_lines[i]=transition_lines[i].replace(s, "")
            transition_lines[i]=transition_lines[i].replace(" ", "")
            
        transition_lines=transition_lines[:-2]
        index1=0
        index2=0
        pop_list=[]
        for l in transition_lines:
            fsa_transitions.append(l)
            if "]]" in l:
                if index1==index2:
                    index1+=1
                    index2+=1
                else:
                    fsa_transitions[index1]=fsa_transitions[index1]+l
                    pop_list.append(transition_lines.index(l))
                    index1=index2+1
                    index2+=1
            elif "[]" in l:
                if index1==index2:
                    index1+=1
                    index2+=1
                else:
                    print("something went wrong at the no transition line")
            else:
                if index1!=index2:
                    fsa_transitions[index1]=fsa_transitions[index1]+l
                    pop_list.append(transition_lines.index(l))
                index2+=1
        if len(pop_list)>0:
            pop_list.reverse()
            for i in pop_list:
                fsa_transitions.pop(i)
        bad_chars=["[", "]"]
        for i in range(len(fsa_transitions)):
            for c in bad_chars:
                fsa_transitions[i]=fsa_transitions[i].replace(c, "")
        for i in range(len(fsa_transitions)):
            fsa_transitions[i]=re.findall("[^,]+,[^,]+", fsa_transitions[i])
        for i in range(len(fsa_transitions)):
            for j in range(len(fsa_transitions[i])):
                temp=fsa_transitions[i][j].split(",")
                temp[0]=int(temp[0])
                temp[1]=int(temp[1])
                fsa_transitions[i][j]=temp
        t=table("sparse", fsa_transitions)
        print(fsa_transitions)
    else:
        t=-1
        print("No formatting assigned to the table object")   
    return t         
def permute_order(fsa, perm):
    """    
        This function permutes the order of the alphabet and the columns
        of the transition table according to the permutation perm, represented
        as a list. The number in position 0 indicates the new index of the 
        first letter of the alphabet, etc.
    """  
    #initialize temp_list with the correct length
    temp_list=[None]*len(perm)
    #reorder alphabet
    for i in range(len(perm)):
        temp_list[perm[i]]=fsa.alphabet.names[i]
    fsa.alphabet.names=temp_list
    #reorder transitions
    for i in range(len(fsa.table.transitions)):
        temp_list=[None]*len(perm)
        for j in range(len(perm)):
            temp_list[perm[j]]=fsa.table.transitions[i][j]
        fsa.table.transitions[i]=temp_list 
        #self.print_fsa("somefilename.wa")    


def compose_auto(fsa1, fsa2):
    """
    Constructs the composed automaton of fsa1 then fsa2. 
    fsa1 and fsa2 should be multiplier 
    automata with the same alphabet since they came from the same general
    multiplier automaton. 
    
    This function works best when fsa1 and fsa2 are 
    deterministic fsa's and have dense tables. 
    We don't need it in another case yet.

    Parameters
    ----------
    fsa1 : dense deterministic fsa object 
    fsa2 : dense deterministic fsa object 

    Returns
    -------
    det fsa object with labeled accept states

    """
    
    #print('fsa1 information')
    #print(fsa1.states.label_names)
    #print(fsa1.states.state_labels)
    #print(fsa1.states.accepting)
    #print('fsa2 information')
    #print(fsa2.states.label_names)
    #print(fsa2.states.state_labels)
    #print(fsa2.states.accepting)
    
    #create new fsa object
    new_alpha=fsa1.alphabet
    #store the tables
    t1=fsa1.table.transitions
    t2=fsa2.table.transitions
    #store the state objects
    s1=fsa1.states
    s2=fsa2.states 
    #product the states
    state_list_base=[]
    accepting_state_base=[]
    #create the lists of pairs of states in a shortlex ordering 
    for i in range(s1.size):
        for j in range(s2.size):
            #In state_list, 
            #the state label is i+1 while the index is i
            state_list_base.append((i+1,j+1))
            if (i+1, j+1) not in accepting_state_base:
                if i+1 in s1.accepting and j+1 in s2.accepting:
                    accepting_state_base.append((i+1, j+1))
    
    #print(state_list_base)
    #print(accepting_state_base)
    
    
    #make new transition table and states
    new_t=[]
    new_state_list=[]
    new_state_trans_done=[]
    #print('making transition table of composed automaton')
    
    #append the start state
    new_state_list.append({(1,1)})
    new_state_trans_done.append(False)
    counter=0
    f = open("compose_multiplier_test.txt", "a")
    
    while False in new_state_trans_done: 
        
        counter=counter+1
        temp_row=[]
        temp_set=set(())
        state_index=new_state_trans_done.index(False)
        curr_state=new_state_list[state_index] #i is a set of tuples {(1,2), (3,4), (5,6)}
        #print(curr_state)
        f.write(str(curr_state)+'\n')
        f.write(str(counter))
        
        #Find the transition function by considering every possible initial
        #tuple i in curr_state 
        #each j in new_alpha.names corresponds to a column of the transition
        #table
        #the union of the outputs over all i's belong to the same row and column
         
        for j in new_alpha.names: #j is a tuple ('a','b')
            f.write("j" + ' ' + str(j)+'\n')
            for i in curr_state:
                f.write("i"+ ' ' + str(i)+'\n')
                #calculate new transition table values when the connecting letter 
                #is not epsilon
                for k in new_alpha.base:
                    fsa1_pair=(j[0],k)
                    fsa2_pair=(k,j[1])
                    index1=new_alpha.names.index(fsa1_pair)
                    index2=new_alpha.names.index(fsa2_pair)
                    
                    #follow along j[0],k from state i[0] in FSA1 
                    #follow along k,j[1] from state i[1] in FSA2
                    #print(fsa1_pair, index1, fsa2_pair, index2, new_tuple)
                    new_tuple=(t1[int(i[0]-1)][index1], t2[int(i[1]-1)][index2])
                    f.write(str(fsa1_pair)+ ', ' +str(index1) + ', '+
                            str(fsa2_pair)+ ', ' +str(index2) + ', '+ 
                            str(new_tuple)+'\n')
                    #only kdeep states that do not contain 0's
                    if new_tuple[0]!=0 and new_tuple[1]!=0:
                        temp_set.add(new_tuple)
                #if there is an epsilon in the tuple, j
                if j[1]!='_' and j[0]!='_':
                    #for padding symbol in 2nd coordinate
                    fsa1_pair=(j[0],'_')
                    index1=new_alpha.names.index(fsa1_pair)
                    fsa2_pair=('_',j[1])
                    index2=new_alpha.names.index(fsa2_pair)
                    new_tuple=(t1[int(i[0]-1)][index1],t2[int(i[1]-1)][index2])
                    #print(fsa1_pair, index1, fsa2_pair, index2, new_tuple)
                    f.write(str(fsa1_pair)+ ', ' +str(index1) + ', '+
                            str(fsa2_pair)+ ', ' +str(index2) + ', '+ 
                            str(new_tuple)+'\n')
                    #only keep states that do not contain 0's
                    if new_tuple[0]!=0 and new_tuple[1]!=0:
                        temp_set.add(new_tuple)  
            #print("temp_set", temp_set)
            f.write("temp_set"+' '+str(temp_set)+'\n')
            
            
            if temp_set != set():
                #append and reset after looping over the tuples in curr_set      
                temp_row.append(temp_set)
                #add temp_set to the list of states
                if temp_set not in new_state_list:
                    new_state_list.append(temp_set)
                    new_state_trans_done.append(False)
            else:
                temp_row.append(0)
            
            temp_set=set()
        #append the row and set the state to "done finding the transition"     
        new_t.append(temp_row)
        new_state_trans_done[state_index]=True
        
        if counter%100==0:
            
            f.write(str(counter)+'\n')
            #print(new_state_list)
            #print(new_state_trans_done)
            #print(new_t)
            num_True=sum(new_state_trans_done)
            total=len(new_state_trans_done)
            f.write(str(len(new_state_list)) + ' ' +
                    str(len(new_state_trans_done)) + ' ' + 
                    str(num_True/total*100) + ' ' +
                    str(len(new_t))+ '\n')
        
    
    f.close()
    #state list should be flattened to consecutive numbers and not the tuples
    #start state is now 1
    #print('new state list', new_state_list)
    f = open('compose_multiplier_states.txt', 'a')
    f.write('new state list'+ ' '+str(new_state_list))
    f.close()
    f = open('compose_multiplier_table.txt', 'a')
    #print(new_t)
    f.write('new transition table'+' '+str(new_t))
    f.close()
    #find new accepting states which have at least one tuple which has
    #an accept state from the input FSA's in each coordinate
    new_accepting_sets=[]
    for state in new_state_list:
        for acc_state in accepting_state_base:
            if acc_state in state and state not in new_accepting_sets:
                new_accepting_sets.append(state)
    #print(new_accepting_sets, "new accepting states")
    #flattened states to be labeld by list_index+1
    #flatten table
    for row in new_t:
        #print("row", row)
        for state in row:
            #print("state", state)
            if state!=0:
                list_index=new_state_list.index(state)
                row_index=row.index(state)
                row[row_index]=list_index+1
    #flatten accepting list
    new_accepting=[]
    for state in new_accepting_sets:
        new_accepting.append(new_state_list.index(state)+1)
    #print(new_accepting)
    
    #pack output and make the fsa
    new_states=states(len(new_state_list), 1, new_accepting)
    new_table=table("dense deterministic", new_t)
    new_fsa=fsa(new_alpha, new_states, new_table)
    return new_fsa

def make_dense(f):
    """changes a sparse deterministic table to a dense deterministic table"""
    sts=f.states
    alpha=f.alphabet
    t=f.table.transitions
    new_t=[]
    if f.table.format=="sparse":
        for st in range(len(t)):
            temp_row=[]
            st_row=str(t[st])
            for j in range(len(alpha.names)):
                #print(st_row, j+1, st, st_row.find("["+str(j+1)+","))
                word="["+str(j+1)+","
                index0=st_row.find(word)
                if index0!=-1:
                    length=len(word)+2 #add two for the comma+space
                    index1=st_row.find("["+str(j+1)+",")
                    sub_st_row=st_row[index1:]
                    index2=sub_st_row.find("]")
                    st_str=sub_st_row[length-1:index2]
                    #print(index1, index2, st_str)
                    int_st=int(st_str)
                else:
                    int_st=""
                #print(j, st, int_st)
                if int_st!="":
                    temp_row.append(int_st)
                else:
                    temp_row.append(0)
            #print(temp_row)
            new_t.append(temp_row)
    else:
        print("table type not supported")
    tbl=table("dense deterministic", new_t)
    temp_fsa=fsa(alpha, sts, tbl)
    temp_fsa.flags=f.flags
    return temp_fsa
    

def make_ndfsa(fsa):
    """change a deterministic fsa into a non-det fsa.
    the fsa in the input may have epsilon transitions that were added in the
    deflation process"""
    sts=fsa.states
    alpha=fsa.alphabet
    t=fsa.table.transitions
    new_t=[]
    if fsa.table.format=="dense deterministic":
    #change each entry of the transition table to a set
        for i in range(len(t)):
            temp=[]
            for j in range(len(t[i])):
                temp.append({t[i][j]})
            new_t.append(temp)
        tbl=ndtable("non-det dense", new_t)
    elif fsa.table.format=="sparse":
        for i in range(len(t)):
            temp=[]
            for j in range(len(t[i])):
                entry=t[i][j]
                temp.append([entry[0], {entry[1]}])
            new_t.append(temp)
        tbl=ndtable("non-det sparse", new_t)
    else:
        print("the table type doesn't work")
    #change state to the set of state in other parts of the fsa
    new_accepting=[]
    for q in sts.accepting:
        new_accepting.append({q})
    sts.accepting=new_accepting
    if len(sts.label_names)>0:
        for i in range(len(sts.state_labels)):
            sts.state_labels[i]=[{sts.state_labels[i][0]},sts.state_labels[i][1]]
    temp_fsa=ndfsa(alpha, sts, tbl)
    temp_fsa.flags.append("non-det dense")
    return temp_fsa

def e_free_ndfsa(nd_fsa):
    """This function takes a non-deterministic fsa with _ transitions, 
    which represent epsilon transitions, and returns a non-deterministic
    fsa without _ transitions. This is implemented for dense tables"""
    alpha=nd_fsa.alphabet
    sts=nd_fsa.states
    tb=nd_fsa.table
    t=nd_fsa.table.transitions
    #make all epsilon transitions one column
    ep_index=[]
    for i in range(len(alpha.names)):
        if (alpha.names[i] == '_'):
            ep_index.append(i)
    if len(ep_index)>1:
        for i in ep_index[1:]:
            for st in range(len(t)):
                t[st][ep_index[0]]=t[st][ep_index[0]].union(t[st][i])
                t[st][ep_index[0]]=rem_0(t[st][ep_index[0]])
    #use the smallest index of epsilon
    epsilon_index=ep_index[0]  
    #delete columns of t for ep_index[1:] and extra alphabet entries
    ep_index.reverse()
    if len(ep_index)>1:
        for i in ep_index[:-1]:
            #print("delete index",i)
            alpha.names.pop(i)
            alpha.size=len(alpha.names)
            for st in range(len(t)):
                t[st].pop(i)        

    #while loop the whole table until taking the e-closure no longer changes it
    new_t=t
    #add the first step of epsilon transitions to e-column
    for st in range(len(t)):
        #add the state for path length 0
        new_set=t[st][epsilon_index].union({int(st+1)})
        new_set=rem_0(new_set)
        #if there is an epsilon transition to another state
        if new_set!=set({int(st+1)}): 
            for q in new_set :
                if q!=0:
                    new_set=new_set.union(t[int(q-1)][epsilon_index])
            new_t[st][epsilon_index]=new_set
            """
            #somehow the reverse edge is being added sometimes.
            elif new_set==set({0}) and x_step(t,t[st][epsilon_index],j,"non-det dense")!=set({0}): #epsilon then letter
                middle_set=t[st][epsilon_index]
                print(st, j, new_set, x_step(t, middle_set, j, "non-det dense"))
                new_set=new_set.union(x_step(t,middle_set,j,"non-det dense"))
                print(new_set)
                new_t[st][j]=new_set
            """
        else:
            new_t[st][epsilon_index]=new_set
    count=1        
    #print("round", count, new_t)
    
    #add all possible epsilon transitions to epsilon column
    while new_t!=e_step(new_t, epsilon_index, len(alpha.names)):
        new_t=e_step(new_t, epsilon_index, len(alpha.names))
        count+=1
        #print("round", count, new_t)
    #remove unneccessary 0's
    for st in range(len(new_t)):
        for j in range(len(alpha.names)):
            new_t[st][j]=rem_0(new_t[st][j])
    #print("e-closure column done", new_t)
    #add the epsilon closure of target states
    new_new_t=[]
    for st in range(len(new_t)):
        new_row=[]
        for j in range(len(alpha.names)-1):
            temp_set=new_t[st][j]
            ep_set=new_t[st][epsilon_index].copy()
            for q in ep_set:
                if t[int(q-1)][j]!={0}:
                    #possible transitions after some number of epsilon transitions
                    temp_set=temp_set.union(t[int(q-1)][j])
            #take epsilon closure of temp_set
            target_ep=set()        
            for q in temp_set:
                if q!=0:
                    #add epsilon closure of states in temp_set 
                    target_ep=target_ep.union(new_t[int(q-1)][epsilon_index])
            temp_set=temp_set.union(target_ep)
            #print(st, j, target_ep, temp_set)
            new_row.append(temp_set)
        new_new_t.append(new_row)            
    #remove unneccessary 0's
    for st in range(len(new_new_t)):
        for j in range(len(alpha.names)-1):
            new_new_t[st][j]=rem_0(new_new_t[st][j])
    
    #finally, delete the epsilon_index in new_t and alphabet
    alpha.names.pop(epsilon_index)
    alpha.size=len(alpha.names)
    tb=ndtable("non-det dense", new_new_t)
    return(ndfsa(alpha, sts, tb))    

def rem_0(q):
    """takes a state and checks if there's an uncessary 0 to remove"""
    if 0 in q:
        q.remove(0)
        if len(q)!=0:
            return q
        else:
            return set({0})
    else:
        return q
    
def e_step(t, ep_index, alpha_len):
    """This function takes the full transition table t (non-det dense)
    and outputs the set of states from taking one epsilon move
    The epsilon step of a state q is the set of states p such that there
    is a single edge from q to p labeled by "_" """
    temp_table=t
    for st in range(len(t)):
        new_set=t[st][ep_index]
        new_set=rem_0(new_set)
        #if there is an epsilon transition to another state
        if new_set!=set({int(st+1)}): 
            for q in new_set :
                if q!=0:
                    new_set=new_set.union(t[int(q-1)][ep_index])
            temp_table[st][ep_index]=new_set
        else:
            temp_table[st][ep_index]=new_set
    return temp_table
 
def x_step(t, st, x_index, table_type):
    """This function takes the full transition table t (non-det dense)
    and outputs the set of states from taking one x move.
    """
    temp_set=set()
    if table_type=="non-det dense":
        for s in st:
            if t[s-1][x_index]!=set({0}):
                temp_set=temp_set.union(t[s-1][x_index])
                #print("x-step", st, s, s-1, x_index, t[s-1][x_index], temp_set)
    else:
        print("table type not supported")
    return temp_set

def make_det(nd_fsa):
    """change a non-det fsa into a det fsa
    Not tested 
    """
    sts=nd_fsa.states
    alpha=nd_fsa.alphabet
    tb=nd_fsa.table
    t=nd_fsa.table.transitions
    new_t=[]
    #apply power set algorithm to transition table
    state_set_list=[{1}] #the initial state is {1} 
    if nd_fsa.table.format=="non-det dense":
        for i in range(len(t)):
            for j in range(len(t[i])):
                if t[i][j] not in state_set_list and t[i][j]!=set({0}): 
                    state_set_list.append(t[i][j])
        #print("state list", state_set_list)
        #create new table
        for st in state_set_list:
            temp_row=[]
            for j in range(len(alpha.names)):
                k=x_step(t,st,j,"non-det dense")
                if k!=set():
                    if k in state_set_list:
                        temp_row.append(state_set_list.index(k)+1)
                    else: 
                        state_set_list.append(k)
                        temp_row.append(state_set_list.index(k)+1)
                else:
                    temp_row.append(0)
            new_t.append(temp_row)
    #TODO test more 
    else:
        print("table type not supported")
    #update table object and change to det table
    tb=table("dense deterministic", new_t)
    #update state object
    sts.size=len(state_set_list)
    
    new_accepting=[]
    for st in state_set_list:
        for q in st:
            #print("temp state ", q)
            if {q} in sts.accepting:
                if state_set_list.index(st)+1 not in new_accepting:
                    new_accepting.append(state_set_list.index(st)+1)
    sts.accepting=new_accepting
    #update initial state to [1]
    init_state = []
    for st in sts.initial[0]:
        init_state.append(st)
    sts.initial = init_state
            
    """
    #I believe this is a carry-over from when I tried to compose automata,
    #but is no longer needed by active code
    if len(sts.label_names)>0: #TODO: figure out what is going on with this
        for i in range(len(sts.state_labels)):
            row=sts.state_labels[i]
            new_state_num=state_set_list.index(row[0])+1
            row[0]=new_state_num
            sts.state_labels[i]=row
    """
    temp_fsa=fsa(alpha, sts, tb)
    temp_fsa.flags=nd_fsa.flags
    return temp_fsa
        
def labeled_min(in_fsa):
    print("Start labeled min function \n\n")
    #find non-accepting states
    nonAccList = []
    AccList = []
    for i in range(in_fsa.states.size):
        if i in in_fsa.states.accepting:
            AccList.append(i)
        else:
            if i!=0:
                nonAccList.append(i)
    nonAccStates = set(nonAccList)
    print("Accepting States", AccList)
    print("Non accepting states", nonAccStates)
    #Find initial state partition
    statePartition = [set({0})]    
    print(in_fsa.states.label_numbers_pairs, '\n\n')
    print(in_fsa.states.label_states_pairs, '\n\n')
    
    for pair in in_fsa.states.label_numbers_pairs:
        statePartition.append(set())
    for pair in in_fsa.states.label_states_pairs:
        statePartition[int(pair[1])].add(int(pair[0]))
    statePartition.append(nonAccStates)
    preRefinedSets = statePartition.copy()
    print("initial state partition", statePartition)
    i=0
    while len(preRefinedSets)!=0:
        #choose and remove a set A from preRefinedSets
        setA = preRefinedSets[0]
        preRefinedSets.pop(0)
        #build X (transitions which reach A)
        X=set()
        for ltr in range(in_fsa.alphabet.size):
            for row in in_fsa.table.transitions:
                if row[ltr] in setA:
                    X.add(row[ltr])       
            #print("letter", ltr, "set reachable to setA", X)
            for setB in statePartition:
                #print(setB.intersection(X), setB.difference(X))
                setBandX = setB.intersection(X)
                setBnotX = setB.difference(X)
                if len(setBandX)!=0 and len(setBnotX)!=0:
                    #print("Inside the if statement")
                    #print("B intersect X", setBandX, "B minus X", setBnotX)
                    if setB in statePartition:
                        indexB = statePartition.index(setB)
                        statePartition.pop(indexB)
                        statePartition.append(setBandX)
                        statePartition.append(setBnotX)
                    else:
                        if len(setBnotX) >= len(setBandX):
                            if setBandX not in preRefinedSets:
                                preRefinedSets.append(setBandX)
                        else:
                            if setBnotX not in preRefinedSets:
                                preRefinedSets.append(setBnotX)
        print("Iterations", i, '\n', "size of partition", len(statePartition), '\n', 
              statePartition)
        i=i+1
    #Create new states object
    newAccList = []
    stateNumPairs = []
    #new accept state list
    for i in range(len(statePartition)):
        for st in statePartition[i]:            
            if st in AccList:
                if i not in newAccList:
                    newAccList.append(i)
                    print("Index", i, "StateSet", statePartition[i])
    #new list of state and label number pairs
    print("find new list of states and label number pairs")
    for i in newAccList:
        stateSet = statePartition[i]
        for pair in in_fsa.states.label_states_pairs:
            if list(stateSet)[0]==int(pair[0]):
                print("state set", stateSet, "current pair", pair, "new pair", [i, pair[1]])
                label = pair[1]
        stateNumPairs.append([i, label])
    #new start state list
    for i in range(len(statePartition)):
        for st in statePartition[i]:
            if st in in_fsa.states.initial:
                newStartState = i
    sts = states(len(statePartition)-1, newStartState, newAccList)
    sts.label_states_pairs = stateNumPairs
    sts.label_numbers_pairs = in_fsa.states.label_numbers_pairs.copy()
    sts.format = 'labeled'
    #Create new transition table
    newTransitions = []
    #statelist without 0 at the top row!
    newStateList = list(range(len(statePartition)))
    newStateList.pop(0)
    print("New state list", newStateList, "start state", newStartState, "accept states", newAccList)
    for i in newStateList:
        row = []
        stateIndex = i-1
        for j in range(in_fsa.alphabet.size):
            oldState = in_fsa.table.transitions[stateIndex][j]
            if oldState == 0:
                newState = 0
            #print("state", i, "alphabet index, j", j, "oldState", oldState)
            else:
                for stateSet in statePartition:
                    if oldState in stateSet:
                        newState = statePartition.index(stateSet)
                    #print("newState", newState)
            row.append(newState)
        newTransitions.append(row)
    tbl = table("dense deterministic", newTransitions)
    print("New Transitions", newTransitions)
    min_fsa = fsa(in_fsa.alphabet, sts, tbl)
    return min_fsa
    
    
#TODO bfs dfsa; is this necessary? 


##############################################################################
# Driver code for testing
"""
fsa_2 = create_fsa_from_file('/home/aurora/Desktop/fsa_unmind')

fsa_1 = labeled_min(fsa_2)

#We see in this example that the minimization algorithm doesn't work.
fsa_3 = create_fsa_from_file('/home/aurora/Desktop/fsa3_unmind')

fsa_4 = labeled_min(fsa_3)
"""

"""
gm_fsa = create_fsa_from_file('/home/aurora/Documents/Diss_Exs/HNN_Exs/gp1/unmingm.gm')
labeled_min(gm_fsa)
"""

"""
a=alphabet(10, ['a','A','b','B','c','C','d','D','t0','T0'], "type?")

tr=[[2,3,4,5,0,0,0,0,6,0],
    [2,0,4,5,0,0,0,0,6,0],
    [0,3,4,5,0,0,0,0,6,0],
    [0,0,4,0,0,0,0,0,6,0],
    [0,0,0,5,0,0,0,0,6,0],
    [0,0,0,0,0,0,7,8,0,0],
    [0,0,0,0,0,0,7,0,0,9],
    [0,0,0,0,0,0,0,8,0,9],
    [0,0,4,5,0,0,0,0,0,0] ]
t=table("dense deterministic", tr)
sts=states(9,[1], [1,2,3,4,5,7,8])
det_fsa=fsa(a, sts, t)
#mimic deflation stuff
for t in ['t0','T0']:
    index_t=det_fsa.alphabet.names.index(t)
    det_fsa.alphabet.names[index_t]="_"
    det_fsa.alphabet.size+=-1
det_fsa.alphabet.size=len(det_fsa.alphabet.names) #to adjust for '_'
#print(det_fsa.alphabet.size)
nd_fsa=make_ndfsa(det_fsa)
print(nd_fsa.alphabet.names)
print(nd_fsa.table.transitions)
efree_nd_fsa=e_free_ndfsa(nd_fsa)
print(efree_nd_fsa.alphabet.names)
print(efree_nd_fsa.table.transitions)
print("try to determinize it")
#defl_fsa=make_det(efree_nd_fsa)
#print(defl_fsa.states.accepting)
#print(defl_fsa.table.transitions)
"""

"""
a=alphabet(4,['0', '1', '2', '_'], 'dense')
tr=[[{1}, {0}, {0}, {2}], [{0}, {2}, {0}, {3}], [{0}, {0}, {3}, {0}]]
t=table("non-det dense", tr)
sts=states(3, [{1}], [{3}])
e_nf=ndfsa(a, sts, t)
e_free_nf=e_free_ndfsa(e_nf)
f=make_det(e_free_nf)
print(f.table.transitions)
print(f.alphabet.names)
print(f.states.accepting)
"""

"""
create_fsa_from_file("/home/aurora/Documents/PythonGroupTheory/Examples/HNN-examples/Z2overZ/v0gp.wa")
f=create_fsa_from_file("/home/aurora/Documents/PythonGroupTheory/Examples/HNN-examples/Z2overZ/v0gp.gm")

create_fsa_from_file("/home/aurora/Documents/PythonGroupTheory/Examples/HNN-examples/Z2overZ/e0gp-f.cos.wa")
create_fsa_from_file("/home/aurora/Documents/PythonGroupTheory/Examples/HNN-examples/Z2overZ/e0gp-f.cos.gm")

composed_auto(f,f)

"""        
#a=alphabet(4, ['a', 'A', 'b', 'B'], "")
#tr=[[2,3,4,5], [2,0,4,5], [0,3,4,5], [0,0,4,0], [0,0,0,5]]
#t=table('dense deterministic', tr)
#sts=states(5, [1], [1,2,3,4,5]) 
#f=fsa(a, sts, t)
#print(f.table.transitions)  
#nf=make_ndfsa(f)
#print(nf.table.transitions)
#det_f=make_det(nf)
#harder test of make_det function
#a=alphabet(4, ['a', 'A', 'b', 'B'], "")
#tr=[[{1,2}, {2,3}, {1,2}, {4}], [{4}, {1,2}, {3,4}, {0}], [{1},{1,2},{1,2},{4}],[{0},{0},{2,3},{4}]]
#t=ndtable('non-det dense', tr)
#sts=states(1, [{1}],[{1}])
#nf_hard=ndfsa(a,sts,t)
#make_det(nf_hard)


#general multiplier for Z2
"""
print("sparse table testing")
a=alphabet(4, ['a', 'A', 'b', 'B'], "product")
if a.type=="product":
    temp_alpha=[]
    for l in a.names:
        for j in a.names:
            temp_alpha.append((l,j))
        temp_alpha.append((l, "_"))
    for l in a.names: 
        temp_alpha.append(("_", l))
    alpha_names=temp_alpha
a.names=alpha_names
a.base=['a', 'A', 'b', 'B']
sparse_tr=[[[1,2],[3,3],[4,4],[5,5],[7,6],[8,7],[9,8],[10,9],
                          [11,10],[12,11],[13,12],[15,13],[16,14],[17,15],
                          [19,16],[20,17],[21,9],[22,5],[23,17],[24,13]],
                         [[1,2],[3,3],[4,4],[5,5],[11,10],[13,12],[15,13],
                          [16,14],[19,16],[20,17],[21,9],[23,17],[24,13]],
                         [[13,3],[15,5]],
                         [[19,4],[20,5]],
                         [],
                         [[7,6],[8,7],[9,8],[10,9],[12,11],[13,12],[15,13],
                          [17,15],[19,16],[20,17],[22,5],[23,17],[24,13]],
                         [[13,7],[15,9]],
                         [[19,8],[20,9]],
                         [],
                         [[13,10],[23,9]],
                         [[13,11],[23,5]],
                         [[13,12],[15,13],[23,17]],
                         [],
                         [[19,14],[24,9]],
                         [[19,15],[24,5]],
                         [[19,16],[20,17],[24,13]],
                         [] 
                        ]
t=table("sparse", sparse_tr)
sts=states(17, [1], [1,2,5,6,9,12,13,16,17])
sts.state_labels=[
                         [1,1],
                         [2,1],
                         [5,2],
                         [6,1],
                         [9,3],
                         [12,1],
                         [13,4],
                         [16,1],
                         [17,5]
                        ]
sts.label_names=[
                    [1,['IdWord']],
                    [2,['A']],
                    [3,['a']],
                    [4,['B']],
                    [5,['b']]
                  ]
gm_vertex=fsa(a,sts,t)
print("starting deterministic table")
#print(gm_vertex.table.transitions)  
den_gm=make_dense(gm_vertex)
#print(den_gm.table.transitions)
test_comp=compose_auto(den_gm,den_gm)
#nf=make_ndfsa(gm_vertex)
#print("non-det table")
#print(nf.table.transitions)
#det_again=make_det(nf)
#compose_auto(gm_vertex,gm_vertex)
"""

"""
#harder sparse
a=alphabet(4, ['a', 'A', 'b', 'B'], "")
tr=[
    [[1,{1,2}],[2,{1}], [3,{3,4}]],
    [[2,{1,2}],[3,{3,4}],[4,{3,4}]],
    [[1,{1}],[2,{1,2}],[3,{1}],[4,{3,4}]] 
    ]
t=ndtable('non-det sparse', tr)
sts=states(1, [1],[1])
nf_hard=ndfsa(a,sts,t)
make_det(nf_hard)
"""
#f.new_letter('t')
#f.permute_order([1,0,2,3,4,5])
#print(f.alphabet.names)
#print(f.table.transitions)
"""
###########################################################################
#newfsa=create_small_fsa(a, 'A')
#print(newfsa.table.transitions)
#newfsa.print_fsa("newfsa.wa")
#######################################################################
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

#test=create_fsa_from_file("Z2_bFirst.cos.wa")
#test.new_letter('t')
#test.print_fsa("Z2_bFirst_full_alpha.cos.wa")
"""

#############################################################################
#Old Code
