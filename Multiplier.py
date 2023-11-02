#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  3 21:45:21 2021

@author: aurora

This file completes the precomputation needed to calculate the general 
multiplier of an HNN extension, and puts it all together according to my 
theorem.
"""

##############################################################################
# imported files
import subprocess
import FSA_manipulations as FSA
import run_auto as run_auto
import Automata_manipulations as Automata
##############################################################################
# global variables
# kbmag_ftn_dir="/home/aurora/gap-4.11.1/pkg/kbmag-1.5.9/standalone/bin/"
############################################################################


class mult_st(object):
    """
    This class is a single state of the multiplier automaton. An instance of 
    a multiplier state is a 5 tuple: (p, q, wd, leftPad, rightPad) where
    p is the state of the left word in the WA of pi_1(G)
    q is the state of the right word in the WA of pi_1(G)
    wd is the word difference found by using the gog.cascade function
    leftPad is 0 if no padding symbol has been read in the first word
            and 1 otherwise
    rightPad is 0 if no padding symbol has been read in the second word
            and 1 otherwise
    """

    def __init__(self, p, q, wd, leftPad, rightPad):
        self.p = p
        self.q = q
        self.wd = wd
        self.leftPad = leftPad #0 or 1
        self.rightPad = rightPad #0 or 1

    def __repr__(self):
        return str((self.p, self.q, self.wd, self.leftPad, self.rightPad))
    

def gm_split(gm, file_dir, kbmag_ftn_dir):
    """
    splits the general multiplier of an automatic structure into the multiplier
    automaton for each generator and epsilon. New fsas are created 

    Parameters
    ----------
    X : general multplier fsa
        its an fsa class with sparse tables, labeled states, and labeled alphabet
    file_dir : working directory

    Returns
    -------
    List of the multiplier FSA's

    """
    print('In GM splitting function')
    X = gm.states.label_numbers_pairs
    Y = gm.states.label_states_pairs
    print("X", X, "type of X[0]", type(X[0]), "type of X[0][0]", type(X[0][0]))
    print("Y", Y)
    multiplier_list = []
    for x in X:
        if type(x[0]) is str: #if read from file it's a string
            label_num = int(x[0].replace(",", ""))
        else: #if created in program, it's an integer
            label_num = x[0]
        label_string = x[1]
        temp_label_names = [[label_num, label_string]]
        print(temp_label_names)
        new_accept = []
        for y in Y:
            if int(y[1]) == label_num:
                new_accept.append(int(y[0]))
        print(label_num, type(label_num), label_string, type(label_string))
        temp_gm = gm
        temp_gm.states.accepting = new_accept
        print("new accepting", new_accept)
        temp_gm.states.label_numbers_pairs = temp_label_names
        temp_gm.states.label_states_pairs = [
            [int(y[0]), int(y[1])] for y in Y if int(y[1]) == label_num]
        print("label_numbers_pairs", temp_gm.states.label_numbers_pairs, 
              "label_states_pairs", temp_gm.states.label_states_pairs)
        # minimize FSA here: note that we can use kbmag's min function since we have the labels stored
        file_name = "ltr"+label_string+".gm"
        temp_gm.print_fsa(file_dir + file_name)
        subprocess.run([kbmag_ftn_dir+"fsalabmin",
                       "-vv", file_dir + file_name])
        # create minimized FSA from file and grab label numbers
        temp_gm = FSA.create_fsa_from_file(file_dir+file_name+".labmin")
        #for i in range(len(temp_gm.states.label_numbers_pairs)):
        #    temp_gm.states.label_numbers_pairs[i][0] = int(
        #        temp_gm.states.label_numbers_pairs[i][0])
        #for i in range(len(temp_gm.states.label_states_pairs)):
        #    temp_gm.states.label_states_pairs[i][0] = int(
        #        temp_gm.states.label_states_pairs[i][0])
        #    temp_gm.states.label_states_pairs[i][1] = int(
        #        temp_gm.states.label_states_pairs[i][1])
        multiplier_list.append([temp_gm, label_string])
    return multiplier_list

def find_all_wds(M, alph):
    "This function finds all words up to length M over the alphabet 'alph'"
    if M >= 0:
        all_wds = set({"IdWord"})
    if M >= 1:
        temp_wds = set()
        for i in alph:
            temp_wds = temp_wds.union({i})
        all_wds = all_wds.union(temp_wds)
    if M >= 2:
        for i in range(2, M+1):
            temp_set = set({})
            for w in all_wds:
                for l in alph:
                    if w != "IdWord" and l != "IdWord":
                        temp_set = temp_set.union({str(w+l)})
            all_wds = all_wds.union(temp_set)
    return all_wds


def wd_check(f, wd):
    """this function checks if the string, wd, is accepted by the deterministic and dense fsa, f"""
    curr_st = 1
    for i in range(len(wd)):
        if wd[i] not in f.alphabet.names:
            return False
        index1 = f.alphabet.names.index(wd[i])
        curr_st = f.table.transitions[int(curr_st-1)][index1]
        if curr_st == 0:
            break
    if curr_st in f.states.accepting:
        return True
    else:
        return False


def pair_wd(M, e, gog, file_dir, kbmag_ftn_dir):
    """This function should return the function  as a pair of ordered lists
    which takes a word,w  and outputs the result, phi_e(w).
    The map phi_e is the isomorphism associated to an edge e of gog.
    """
    # Find all words up to particular length M
    alph_f = e.fg.gens+e.fg.invgens
    # find all words in the coset auto structure for the pair G_v and G_e(forward)
    all_f_wd = find_all_wds(M, alph_f)

    # check if each word is in the vertex wa
    tv_wa = gog.vertices[e.tv].group.wa
    good_f_wd = set({})
    for wd in all_f_wd:
        if wd == "IdWord":
            good_f_wd.add(wd)
        else:
            if wd_check(tv_wa, wd) == True:
                good_f_wd.add(wd)
    list_good_f_wd = []
    list_good_r_wd = []
    for i in good_f_wd:
        list_good_f_wd.append(i)
        if i != "IdWord":
            list_good_r_wd.append(i.translate(e.isom))
        else:
            list_good_r_wd.append("IdWord")

    print("list of good forward words", list_good_f_wd)
    print("list of good reverse words", list_good_r_wd)
    wd_ftn = (list_good_f_wd, list_good_r_wd)
    return wd_ftn


def new_pair_wd(mid_wd, st_ltr, pair_wd_ftns):
    """
    From the middle word and the stable letter, this function finds the new 
    paired word which is calculated from the isomorphism provided by the user. 

    """

    print("In new paired word function")
    print(pair_wd_ftns[0][0])
    print(st_ltr)
    print(mid_wd)
    if st_ltr == 't0':
        if mid_wd in pair_wd_ftns[0][0]:
            index1 = pair_wd_ftns[0][0].index(mid_wd)
            print(pair_wd_ftns[0][1][index1])
            new_mid_wd = pair_wd_ftns[0][1][index1]
        else:
            return 0
    else:
        # st_ltr='T0'
        if mid_wd in pair_wd_ftns[0][1]:
            index1 = pair_wd_ftns[0][1].index(mid_wd)
            print(pair_wd_ftns[0][0][index1])
            new_mid_wd = pair_wd_ftns[0][0][index1]
        else:
            return 0
    return new_mid_wd


def bdd_wd_machs(M, gog, obj, obj_type, file_dir, kbmag_ftn_dir):
    """
    Computes the bounded word machines for all accepted words up to length M
    for the vertex group, for the forward edge subgroup's coset auto structure 
    and for the reverse edge subgroup's coset auto structure. 

    Parameters
    ----------
    M : integer
    max bound on length of words
    obj : either a vertex or an edge class
    obj_type: 'vertex' or 'edge'
    file_dir : string
        the absolute file directory where the FSA files and GOG files exist
    kbmag_ftn_dir : string
        the absolute file directory where the kbmag binary function files exist

    Returns
    -------
    If obj_type is vertex, then a list of lists of lists is returned.  
    The items of the first layer are [list of word multipliers, vX] where
    X is the label of the vertes. 
    Each item of the list of word multipliers is the list 
    [multiplier automata, word difference] where the multiplier automata
    accepts words in normal form with the specified word difference.

    If obj_type is edge, then a list of lists of lists is returned. 
    The items of the first layer are 
    The first entry is the forward group's coset representatives of length at most M, 
    the second entry is the bounded word machines (fsa's) of the words in the first entry,
    the third entry is the reverse group's coset representatives of length at most M,
    and the fourth entry is the bounded word machines (fsa's) of the words in the third entry.
    """
    print("In bounded word machine function for", obj_type)
    # find all words of length less than M and convert to a list

    if obj_type == 'vertex':
        list_good_wds = []
        all_wds_v0 = find_all_wds(M, obj.group.gens)
        v0_wa = obj.group.wa
        good_wds_v0 = set({})
        for wd in all_wds_v0:
            if wd == 'IdWord':
                good_wds_v0.add(wd)
            else:
                if wd_check(v0_wa, wd) == True:
                    good_wds_v0.add(wd)
        for wd in good_wds_v0:
            list_good_wds.append(wd)
        print(list_good_wds)
    elif obj_type == 'edge':
        # find the alphabet of the coset structures
        ef_coswa = obj.fg.coswa
        tv_alpha = gog.vertices[obj.tv].group.gens + \
            gog.vertices[obj.tv].group.invgens
        ef_coswa_alpha = [a for a in tv_alpha if a not in obj.fg.gens and
                          a not in obj.fg.invgens]
        er_coswa = obj.rg.coswa
        iv_alpha = gog.vertices[obj.iv].group.gens + \
            gog.vertices[obj.iv].group.invgens
        er_coswa_alpha = [a for a in iv_alpha if a not in obj.rg.gens and
                          a not in obj.rg.invgens]
        # find all words
        all_f_wds_e = find_all_wds(M, ef_coswa_alpha)
        all_r_wds_e = find_all_wds(M, er_coswa_alpha)

        # print(all_f_wds_e)
        # print(all_r_wds_e)

        # find all good words by using the coset word acceptor
        good_wds_ef = set({})
        for wd in all_f_wds_e:
            if wd == 'IdWord':
                good_wds_ef.add(wd)
            else:
                if wd_check(ef_coswa, wd) == True:
                    good_wds_ef.add(wd)
        good_wds_er = set({})
        for wd in all_r_wds_e:
            if wd == 'IdWord':
                good_wds_er.add(wd)
            else:
                if wd_check(er_coswa, wd) == True:
                    good_wds_er.add(wd)

        # turn the sets into lists
        list_good_f_wds = []
        for wd in good_wds_ef:
            list_good_f_wds.append(wd)
        list_good_r_wds = []
        for wd in good_wds_er:
            list_good_r_wds.append(wd)

    else:
        return "Error, you did not provide a vertex or edge flag for the object."
    # HNN specific coding
    # isolate the list of multipliers with labels
    if obj_type == 'vertex':
        v0_mult_list = obj.group.mult_list
        # print(v0_mult_list)
        # print(list_good_wds)
    elif obj_type == 'edge':
        ef_mult_list = obj.fg.mult_list
        er_mult_list = obj.rg.mult_list

    # An entry of word_mult_list is a list whose first entry is the string, x
    # the second entry is the automata accepting pairs of words (v,w) satisfying
    # vx=w in the group or coset structure

    """
    #Testing code:
    temp_auto=FSA.compose_auto(v0_mult_list[2][0], v0_mult_list[2][0])
    """

    if obj_type == 'vertex':
        word_mult_list = []
        for wd in list_good_wds:
            # convert the word as a string into a list
            if wd != "IdWord":
                wd_list = list(wd)
            else:
                wd_list = ['IdWord']
            # obtain the labels of the list of multipliers
            label_slice = [a[1] for a in v0_mult_list]
            idword_index = label_slice.index('IdWord')
            # make the starting automaton M_{IdWord}
            temp_auto = v0_mult_list[idword_index][0]
            if wd != 'IdWord':
                # iteratively compose the automaton one letter at a time to form
                # the word multiplier
                for i in range(len(wd_list)):
                    index = label_slice.index(wd_list[i])
                    #print('in for loop,', v0_mult_list[index][1])
                    temp_auto = FSA.compose_auto(
                        temp_auto, v0_mult_list[index][0])
            word_mult_list.append([temp_auto, wd])
        # print(word_mult_list)
    elif obj_type == 'edge':
        f_word_mult_list = []
        r_word_mult_list = []
        for wd in list_good_f_wds:
            # convert the word as a string into a list
            if wd != "IdWord":
                wd_list = list(wd)
            else:
                wd_list = ['IdWord']
            # obtain the labels of the list of multipliers
            f_label_slice = [a[1] for a in ef_mult_list]
            for label in f_label_slice:
                if "IdWord" in label:
                    idword_index = f_label_slice.index(label)
            # make the starting automaton M_{IdWord}
            temp_auto = ef_mult_list[idword_index][0]
            if wd != 'IdWord':
                # iteratively compose the automaton one letter at a time to form
                # the word multiplier
                for i in range(len(wd_list)):
                    index = f_label_slice.index(wd_list[i])
                    temp_auto = FSA.compose_auto(
                        temp_auto, ef_mult_list[index][0])
            f_word_mult_list.append([temp_auto, wd])
        for wd in list_good_r_wds:
            # convert the word as a string into a list
            if wd != "IdWord":
                wd_list = list(wd)
            else:
                wd_list = ['IdWord']
            # obtain the labels of the list of multipliers
            r_label_slice = [a[1] for a in er_mult_list]
            for label in r_label_slice:
                if "IdWord" in label:
                    idword_index = r_label_slice.index(label)
            # make the starting automaton M_{IdWord}
            temp_auto = er_mult_list[idword_index][0]
            if wd != 'IdWord':
                # iteratively compose the automaton one letter at a time to form
                # the word multiplier
                for i in range(len(wd_list)):
                    index = r_label_slice.index(wd_list[i])
                    temp_auto = FSA.compose_auto(
                        temp_auto, er_mult_list[index][0])
            r_word_mult_list.append([temp_auto, wd])

    # output
    if obj_type == 'vertex':
        return [word_mult_list, 'v0']
    elif obj_type == 'edge':
        return [f_word_mult_list, 'fe'+str(obj.label),
                r_word_mult_list, 're'+str(obj.label)]

# new bdd word


def new_bdd_wd(mid_wd, ltr_pair, vfile, kbmag_ftn_dir, file_dir):
    # TODO rewrite this function to call vx.group.word_reduce function
    # TODO move to coding graveyard 
    # The new middle word defined by using word reduce in the appropriate
    # vertex group's automatic structure
    l = ltr_pair[0]
    r = ltr_pair[1]
    new_word = str()
    if l != "IdWord":
        new_word = new_word+l.swapcase()
    if mid_wd != "IdWord":
        for i in range(len(mid_wd)):
            if new_word != str():
                new_word = new_word+"*"+mid_wd[i]
            else:
                new_word = mid_wd[i]
    if r != "IdWord":
        if new_word != str():
            new_word = new_word+"*"+r
        else:
            new_word = r
    if new_word == str():
        new_word = "IdWord"
    new_word = new_word+";"
    output = subprocess.run([kbmag_ftn_dir+"wordreduce", file_dir + vfile],
                            input=new_word.encode(), capture_output=True)
    print("inputs:", "mid_wd =", mid_wd, "ltr_pair = (",
          l, ",", r, ")", "as product: ", new_word)
    print("output:")
    txt_out = str(output.stdout)
    if "reduces to:" in txt_out:
        index1 = txt_out.index("reduces to:")
        suffix = txt_out[index1+15:]
        index2 = suffix.index("\\")
        suffix = suffix[:index2]
        print("suffix", suffix, "type", type(suffix))
        new_mid_wd = str()
        if suffix == "IdWord":
            return "IdWord"
        while suffix != str():
            if "*" in suffix:
                index3 = suffix.index("*")
                temp = suffix[:index3]
                print(suffix, "   ", temp)
                if "^" in temp:
                    index4 = temp.index("^")
                    letter = temp[:index4]
                    print(letter)
                    number = int(temp[index4+1:])
                    for _ in range(number):
                        new_mid_wd = new_mid_wd+letter
                else:
                    new_mid_wd = new_mid_wd+temp
                suffix = suffix[index3+1:]
                print("new suffix", suffix)
            else:
                temp = suffix
                if "^" in suffix:
                    index4 = temp.index("^")
                    letter = temp[:index4]
                    print(letter)
                    number = int(temp[index4+1:])
                    for _ in range(number):
                        new_mid_wd = new_mid_wd+letter
                    suffix = str()
                else:
                    new_mid_wd = new_mid_wd+suffix
                    suffix = str()
    else:
        print(txt_out)
        new_mid_wd = 'ERROR: word reduce function failed'
    print("new middle word", new_mid_wd)
    return new_mid_wd

def testMultiplier(testFSA, fsaString, goalWAString, file_dir, kbmag_ftn_dir, out_dir):
    #This function tests if the projection of the multiplier for fsaString 
    #matches the goalWA
    print("Test multiplier for ", fsaString)
    testFSA.print_fsa("testingTestMultiplier" + fsaString + ".txt")
    eps_alpha = testFSA.alphabet.base.copy()
    eps_alpha.append("_")
    letter_alphabet = FSA.alphabet(len(eps_alpha), eps_alpha, "identifiers")
    #The intial state {1} is always an accept state since we are
    #trying to accept the pairs (_, x) for x in X
    test_acc_states = [{1}]

    test_table=[]
    isAcceptState = True
    #find the transition table for the projection onto the first coordinate
    for iRow in range(len(testFSA.table.transitions)):
        temp_row = []
        for iFirstLetter in range(len(testFSA.alphabet.base)+1):
            new_st = set()
            for iSecondLetter in range(len(testFSA.alphabet.base)+1):
                alpha_index = iFirstLetter*(len(testFSA.alphabet.base)+1)+iSecondLetter
                #print("letter index testing.")
                #print("first letter: ", iFirstLetter, " second letter: ", iSecondLetter)
                #print("alphabet index: ", alpha_index)
                if alpha_index != len(eps_alpha)*len(eps_alpha)-1:
                    new_st.add(testFSA.table.transitions[iRow][alpha_index])
            new_st = FSA.rem_0(new_st)
            
            #check if it's an accept state
            for iState in new_st:
                if iState not in testFSA.states.accepting:
                    isAcceptState = False
            if new_st != {0}:
                print("state: ", new_st, "isAccept: ", isAcceptState)
                print("check if statement1: ", new_st not in test_acc_states)
            if isAcceptState == True:
                #check if in accept state list
                if new_st not in test_acc_states:
                    test_acc_states.append(new_st)
            #reset to True
            isAcceptState = True
            temp_row.append(new_st)
        test_table.append(temp_row)
    print("test table info", len(test_table))
    print("new test table", test_table)   
    print("accept states", test_acc_states)
    #define state object
    letter_states = FSA.states(len(test_table), [{1}], test_acc_states)
    #define table object
    letter_table_object = FSA.table("non-det dense", test_table)
    #define the ndfsa with epsilon transitions
    letter_ndfsa_with_eps = FSA.ndfsa(letter_alphabet, letter_states, letter_table_object)
    #make it an ndfsa without epsilon transitions
    letter_ndfsa = FSA.e_free_ndfsa(letter_ndfsa_with_eps)
    #make it a det fsa 
    letter_fsa = FSA.make_det(letter_ndfsa)
    #print out fsa file and compare the projection to the NF
    projection_file = fsaString+".txt"
    letter_fsa.print_fsa(file_dir + projection_file)
    output = subprocess.run([kbmag_ftn_dir+"fsalequal", out_dir + goalWAString, file_dir+projection_file])
    #print("fsaequal output for letter: "+ fsaString + " " + str(output.returncode)+ "\n")
    return output.returncode

def find_k(gog, travel_info, file_dir, kbmag_ftn_dir, verbose):
    Manual = False
    if Manual == True:
        # construct the fellow traveling constant. Call it M
        if len(travel_info)==0:
            M = find_M(gog, file_dir, kbmag_ftn_dir)
        elif len(travel_info)==1:
            M = int(travel_info[0])
        else:
            #calculating n1
            Ge_wa = FSA.create_fsa_from_file(file_dir + gog.vertices[0].subgp+".wa")
            Ge_gm = FSA.create_fsa_from_file(file_dir + gog.vertices[0].subgp+".gm")
            Ge_gm_list = gm_split(Ge_gm, file_dir, kbmag_ftn_dir)
            gm_sizes = []
            for gm in Ge_gm_list:
                gm_sizes.append(gm[0].states.size)
            n_gm=max(gm_sizes)
            print(gm_sizes)
            n_ge=Ge_wa.states.size
            print(n_ge)
            n1 = max(n_ge, n_gm)+ 1
            #grab fellow travel info
            k2 = int(travel_info[0])
            K = int(travel_info[1])

            k3=K^2 + 1
            print('k2', k2, 'K', K, 'n1', n1, 'k3', k3)
            #calculating M
            k4 = max(k3, k2)
            M = k4*k4 + k4*n1
        print("Fellow traveling constant", M)
    return M

def make_gm(gog, travel_info, file_dir, kbmag_ftn_dir, out_dir, verbose):
    # vertex group general multiplier and edge group coset general multipliers
    # We need the aut struct of v0 and the coset aut structure for every
    # vertex group and edge subgroup pair

    # Note only for HNN extensions
    # v now represents v0
    v = gog.vertices[0]

    # read in gm files
    v.group.gm = FSA.create_fsa_from_file(file_dir + v.group.gm_file)
    v.group.gm.flags.append("gm")
    for e in gog.edges:
        e.fg.cosgm = FSA.create_fsa_from_file(file_dir + e.fg.cosgm_file)
        e.rg.cosgm = FSA.create_fsa_from_file(file_dir + e.rg.cosgm_file)
        e.fg.cosgm.flags.append("cosgm")
        e.rg.cosgm.flags.append("cosgm")

    # create the alphabet of gm
    gm_alph = FSA.alphabet(0, "", "product")
    gm_alph.base = gog.wa.alphabet.names
    X = []
    for x in gm_alph.base:
        for y in gm_alph.base:
            X.append(tuple((x, y)))
        X.append(tuple((x, "_")))
    for x in gm_alph.base:
        X.append(tuple(("_", x)))
    gm_alph.size = len(X)
    gm_alph.names = X

    print("**************************************************")
    print("general multiplier alphabet")
    print("**************************************************", "\n")
    #for each letter of the word acceptor alphabet and Epsilon
    #find the k that makes the multiplier the same as the normal
    #form language
    #In this process, save the transition table and states of 
    #each gm that is created
    
    letter_TF_list = [False] #initial false for identity
    letter_fellowConstant = [0] #initial constant for identity 
    for iLetter in gm_alph.base:
        letter_TF_list.append(False) #false for each letter  
        letter_fellowConstant.append(0) #0 for each letter
    generalMultiplierList = []
    multiplierLetters = ["IdWord"] + gm_alph.base.copy() 
    
    k=0
    #start the while loop for increasing k
    #while k < 4:
    while False in letter_TF_list:
        tempFSA = draft_multiplier(k, gm_alph, gog, file_dir, kbmag_ftn_dir, out_dir, verbose)
        tempFSA.print_fsa(file_dir + "gmFellowConstant"+str(k)+".txt")
        generalMultiplierList.append(tempFSA)
        for index in range(len(letter_TF_list)):
            if letter_TF_list[index] == False:
                if index < len(tempFSA.states.label_numbers_pairs):
                    acceptingForLetter = []
                    for iState in range(len(tempFSA.states.label_states_pairs)):
                        if tempFSA.states.label_states_pairs[iState][1]==tempFSA.states.label_numbers_pairs[index][0]:
                            acceptingForLetter.append(tempFSA.states.label_states_pairs[iState][0])
                    tempFSA = FSA.fsa(generalMultiplierList[k].alphabet, generalMultiplierList[k].states, generalMultiplierList[k].table)
                    tempFSA.states.accepting = acceptingForLetter
                    print("Test for the letter ", multiplierLetters[index])
                    print("temp FSA accepting", tempFSA.states.accepting)
                    print("general multiplier accepting", generalMultiplierList[k].states.accepting)
                    
                    testResult = testMultiplier(tempFSA, multiplierLetters[index], gog.wa_file, file_dir, kbmag_ftn_dir, out_dir)
                    print("Test result for the letter: ", multiplierLetters[index], testResult)
                    if testResult == 0:
                        letter_TF_list[index] = True
                        letter_fellowConstant[index] = k
                
        k+=1
        print("Increasing k to ", str(k))
        print("Current letter T/F list: ", letter_TF_list)
    
    #find the maximum of letter_fellowConstant
    calculatedK = max(letter_fellowConstant)
    print("The fellow traveling constant that works for every letter is ", calculatedK)
    #Finally, create unminimized gm object
    unmin_gm = generalMultiplierList[calculatedK]
    unmin_gm.flags.append("gm")
    unmin_gm.print_fsa(file_dir+"unmingm.gm")
    gm_states = unmin_gm.states
    #For each accepting state label, Minimize unmin_gm using the Automata package
    tableTranspose = Automata.table_transpose(unmin_gm.table.transitions, unmin_gm.alphabet.size, unmin_gm.states.size)
    for iLetter in range(len(gm_states.label_numbers_pairs)):
        acceptingForLetter = []
        for iState in range(len(gm_states.label_states_pairs)):
            if gm_states.label_states_pairs[iState][1]==gm_states.label_numbers_pairs[iLetter][0]:
                acceptingForLetter.append(gm_states.label_states_pairs[iState][0])
            if verbose == True:
                print("Minimizing piece of GM which accepts (w,v) for the letter "+ str(gm_states.label_numbers_pairs[iLetter][1]))
                print(acceptingForLetter)
        gm_aut = Automata.automata("det", 
                                unmin_gm.states.size, 
                                unmin_gm.alphabet.size, 
                                unmin_gm.states.initial, 
                                acceptingForLetter, 
                                tableTranspose)
        file_name = "MinimizeAutomataScriptLetter_"+str(gm_states.label_numbers_pairs[iLetter][1])+".g"
        output_file_name = "autminLetter_"+str(gm_states.label_numbers_pairs[iLetter][1])+".txt"
        Automata.gap_script(gm_aut, file_dir, file_name, output_file_name)
        subprocess.run(["gap", file_dir+file_name])
        numStates, initState, acceptStates, fsaTable = Automata.read_automata_from_file(file_dir, output_file_name)
        if verbose == True:
            print("Number of States", numStates)
            print("Initial State", initState)
            print("Accept States", acceptStates)
            print("FSA Table", fsaTable)
        
        
        #Create state object for each multiplier
        mul_states=FSA.states(numStates, initState, acceptStates)
        mul_states.format="labeled"
        mul_states.label_numbers_pairs= gm_states.label_numbers_pairs
        mul_states.label_states_pairs= [[st, int(iLetter+1)] for st in acceptStates]
        if verbose == True:
            print("label and number pairings", mul_states.label_numbers_pairs)
            print("number and state pairings", mul_states.label_states_pairs)
        #Create table object for each multiplier
        mul_table = FSA.table("dense deterministic", fsaTable)     
        #Create unminimized gm
        tempFSA = FSA.fsa(gm_alph, mul_states, mul_table)
        gog.mult_list.append([str(gm_states.label_numbers_pairs[iLetter][1]), tempFSA] )
    
    # full_gm = FSA.labeled_min(unmin_gm)
    # full_gm.print_fsa(file_dir + 'HigginsNF.gm')
    return 0

def draft_multiplier(k, gm_alph, gog, file_dir, kbmag_ftn_dir, out_dir, verbose):
    
    #create fail state
    fail_state = mult_st(0, 0, 0, 0, 0)
    state_list = [fail_state]
    state_done = [True]
    # create start state
    start_state = mult_st(1, 1, "IdWord", 0, 0)
    state_list.append(start_state)
    state_done.append(False)
    acc_states = [state_list.index(start_state)]

    # append entries of the transition table to transition_test
    print("**************************************************")
    print("transition function process for k set to "+str(k))
    print("**************************************************")
    if verbose ==True:
        f = open(file_dir+"transition_test.txt", "w")
    # build transition function and state list
    # only add accessible states to the list of states
    new_t = []
    IsAcceptState = False
    while False in state_done:
        new_row = []
        index = state_done.index(False)
        s = state_list[index]
        for (x, y) in gm_alph.names:
            new_s , IsAcceptState = find_trans(gog, s, (x, y), k, kbmag_ftn_dir, file_dir)
        
            if verbose == True:
                f.write("current state ")
                f.write("("+str(s.p) + ", "+str(s.q)+", " + str(s.wd)+", "+str(s.leftPad)+", "+str(s.rightPad)+"), ")
                f.write("letter pair ")
                f.write("("+str(x)+","+str(y)+")"+", ")
                f.write("next state ")
                f.write("("+str(new_s.p)+", "+str(new_s.q)+", "+str(new_s.wd)+", "+str(new_s.leftPad)+", "+str(new_s.rightPad)+"), ")
            # check if new_s is in the state_list
            in_list = False
            state_index = "TBA"
            for state in state_list:
                if int(new_s.p) == int(state.p) and int(new_s.q) == int(state.q) and new_s.wd == state.wd and int(new_s.leftPad)==int(state.leftPad) and int(new_s.rightPad)==int(state.rightPad):
                    in_list = True
                    state_index = state_list.index(state)
        #print(str(new_s), "in list", in_list, "with index ", state_index)
        #print("types and length of new state:", type(new_s.p), type(new_s.q), type(new_s.wd))
        
        #Otherwise, add new_s to the state list unless it is the fail state        
            if in_list == False and new_s.wd != 0:  # do not add the fail state to the state list
                state_list.append(new_s)
                state_done.append(False)
                state_index = len(state_list)-1
                #if IsAcceptState is True, add index(new_s) to the accepting state list
                if IsAcceptState == True:
                    acc_states.append(state_index)
            
            if verbose == True:
                #print("index of new state", state_index)
                #print("The state list: ", str(state_list))
                f.write("index in state_list ")
                f.write(str(state_index))
                f.write("Accepting: " + str(IsAcceptState)+"\t" )

            #the state is identified by the index in state_list
            if new_s.wd != 0:
                state_label = int(state_index)
            else: #the label of the fail state is 0
                state_label = int(0)
            #add the entry to the current row of the transition table
            new_row.append(state_label)
            if verbose == True:
                #print("state list: ", str(state_list))
                #print("\n\n\n")
                f.write("new state:"+ str(new_s)+ " new state label: "+ str(state_label)+ "\n")
        #add the current row to the new transition table
        new_t.append(new_row)
        if verbose == True:
            f.write("state list: "+ str(state_list)+ "\n")
            f.write(str(new_row))
            f.write("\n")
        state_done[index] = True
        #end of while for finding states
    if verbose == True:
        f.close()

    #Save the labels of the states and find unique labels
    label_list = [state_list[i].wd for i in range(len(state_list)) ]
    acc_wd_list = [state_list[int(acc_states[i])].wd for i in range(len(acc_states))]
    unique_acc_labels = []
    for i in range(len(acc_wd_list)):
        if acc_wd_list[i] not in unique_acc_labels:
            unique_acc_labels.append(acc_wd_list[i])
    label_numbers = [[i+1, unique_acc_labels[i]] for i in range(len(unique_acc_labels))]
    label_state_pairs = [[acc_states[i], unique_acc_labels.index(acc_wd_list[i])+1] for i in range(len(acc_states))]
    if verbose == True:
        print("state list has length ", str(len(state_list)))
        print("state list", state_list)
        print("number of accepting states", str(len(acc_states)))
        print("accept state list", acc_states)
        print("accept word differences", acc_wd_list)
        print("number of rows of new transition table", str(len(state_list)))
        print("list of word differences", label_list, "with length", len(label_list))
        print("label state pairs for printing", label_state_pairs)
        print("label numbering", label_numbers)
    
    #Convert state labels to index in state_list but pop 0
    state_list = [int(i) for i in range(len(state_list)) ]
    state_list.pop(0)

    if verbose == True:
        print(state_list) 
        print("state list length after discarding fail state", len(state_list))

    #Make GM state labels and state object
    gm_states = FSA.states(len(state_list), [state_list[0]], acc_states)
    gm_states.format="labeled"
    gm_states.label_numbers_pairs = label_numbers
    gm_states.label_states_pairs = label_state_pairs
    #Create table object for gm
    gm_table = FSA.table("dense deterministic", new_t)
    testFSA = FSA.fsa(gm_alph, gm_states, gm_table)
    return testFSA

def find_trans(gog, s, ltr_pair, M, kbmag_ftn_dir, file_dir):
    """
    given a multiplier state s and an element of the alphabet (x,y) find the 
    appropriate transition based on my dissertation!

    Input:
        gog is the graph of groups

        s is a state comprised of the 5-tuple (p, q, wd, padLeft, padRight) where
            p is a state of wa
            q is a state of wa
            wd is a word over the GoG alphabet representing the word difference
            padLeft is a boolean 0: no, pad symbol on the left;
                                1: yes, pad symbol on the left;
            padRight is a boolean 0: no, pad symbol on the right;
                                1: yes, pad symbol on the right;
                                
            OLD: stab is the last seen stable letter ???Is this needed?


        ltr_pair is a list of strings and it should have length 2. 

        M is the integer which determines the maximum length of a word difference
            for new_s.

    Returns the new state, new_s
    """
    wa = gog.wa
    a_names = wa.alphabet.names
    tbl = wa.table.transitions
    wa_acc_states = wa.states.accepting
    
    IsAcceptState = False
    # find new p and q
    #pairs (x,y)
    if ltr_pair[0] != "_" and ltr_pair[1] != "_":
        if s.rightPad==1 or s.leftPad==1:
            return mult_st(0,0,0,0,0), False
        index1 = a_names.index(ltr_pair[0])
        index2 = a_names.index(ltr_pair[1])
        new_s = mult_st(int(tbl[int(s.p-1)][index1]), int(tbl[int(s.q-1)][index2]), "IdWord", s.leftPad, s.rightPad)
        #print(" letter pair ", ltr_pair," indices of letter pair ", str((index1, index2)), " current state ", s, "new state ", new_s)
    #pairs (_,y)
    elif ltr_pair[0] == "_":
        if s.rightPad==1:
            return mult_st(0,0,0,0,0), False
        index2 = a_names.index(ltr_pair[1])
        new_s = mult_st(s.p, tbl[int(s.q-1)][index2], "IdWord", 1, s.rightPad)
    #pairs (x,_)
    else:
        if s.leftPad==1:
            return mult_st(0,0,0,0,0), False
        index1 = a_names.index(ltr_pair[0])
        new_s = mult_st(tbl[int(s.p-1)][index1], s.q, "IdWord", s.leftPad, 1)
    # check for the fail state in either coordinate
    #print("Checking new state tuple \n", "first coordinate", new_s.p, "is it =0?", new_s.p ==0, 
    #      "\n second coordinate", new_s.q, "it it =0?", new_s.q == 0)
    if new_s.p == 0 or new_s.q == 0:
        return mult_st(0, 0, 0, 0, 0), False
    #Check for accept state depending on WA
    if new_s.p in wa_acc_states and new_s.q in wa_acc_states:
        IsAcceptState = True
    

    # find new wd
    # (x,_)
    if "_" not in ltr_pair[0] and "_" in ltr_pair[1]:
        if s.wd == "IdWord":
            new_s.wd = casc_fred_combo(ltr_pair[0].swapcase(), gog,
                                   file_dir, kbmag_ftn_dir)
        else:
            new_s.wd = casc_fred_combo(ltr_pair[0].swapcase()+s.wd, gog,
                                   file_dir, kbmag_ftn_dir)
    # (_,y)
    elif "_" in ltr_pair[0] and "_" not in ltr_pair[1]:
        if s.wd == "IdWord":
            new_s.wd = casc_fred_combo(ltr_pair[1], gog,
                                   file_dir, kbmag_ftn_dir)
        else:
            new_s.wd = casc_fred_combo(s.wd + ltr_pair[1], gog,
                                   file_dir, kbmag_ftn_dir)
    # (x,y)
    else:
        if s.wd == "IdWord":
            new_s.wd = casc_fred_combo(ltr_pair[0].swapcase() + ltr_pair[1], gog,
                                   file_dir, kbmag_ftn_dir)
        else:
            new_s.wd = casc_fred_combo(ltr_pair[0].swapcase()+s.wd+ltr_pair[1], gog,
                                   file_dir, kbmag_ftn_dir)
       # if ltr_pair[0] == 'a' and ltr_pair[1] == 'A':
       #     print(new_s)
    #Check for length for accept state versus standard state
    if wd_len(new_s.wd) >1:
        IsAcceptState = False 
    #Check for fail state
    if wd_len(new_s.wd) > M:
        IsAcceptState = False
        new_s = mult_st(0, 0, 0, 0, 0)

    return new_s, IsAcceptState


def wd_len(word):
    """
    Determine the length of a word by deleting the numbers and calculating
    the length of the new string.

    Parameters
    ----------
    word : string over the alphabet (A_v \cup {stable letters})^{+/- 1}

    Returns
    -------
    n : the length of the alphabet (A_v \cup {stable letters})^{+/- 1}

    """
    temp = word
    nums = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
    #find length of string without the numbers 
    for i in nums:
        temp = temp.replace(i, str())
    n = len(temp)
    #check for IdWord
    if word =="IdWord":
        n=0
    return n

def casc_fred_combo(w, gog, file_dir, kbmag_ftn_dir):
    """
    This function combines the cascading and free reduction functions. 
    It calls a free reduction at the end of each cascade 
    process and if there was a free reduction, the cascade has to be 
    completed again.

    Parameters
    ----------
    w: string representing a word in the fundamental group of the gog
    
    gog : Graph of groups class
    
    file_dir : full directory of the working folder
    
    kbmag_ftn_dir : full directory of the kbmag bin folder

    Returns
    -------
    a word as a string which is in normal form according to the cascade process
    and has letter*(letter inverse) pairs appropriately reduces to the empty
    string or IdWord.

    """
    IsFreeReduced = False
    IsCascaded = False
    temp_word = w
    counter = 0
    
    while IsFreeReduced == False or IsCascaded == False:
        #print("The counter for the while loop is ", str(counter))
        if IsFreeReduced == False:
            #if w is not free reduced, perform that function
            #print("Perform a free reduction on", temp_word, "which has length", len(temp_word))
            new_temp_word = gog.wa.alphabet.free_reduction(temp_word)
            #print("after free reduction, the temp word is ", new_temp_word, "which has length", len(temp_word))
            IsFreeReduced = True
            #Special Case: if only a free reduction needs to happen#
            #if temp_word == str():
            #    return "IdWord"
            #if the word changed, cascading needs to be performed again
            #and temp_word needs to be updated
            if new_temp_word != temp_word:
                IsCascaded = False
                temp_word = new_temp_word
        if IsCascaded == False:
            #if w is not cascaded, perform that function
            #print("Perform a cascade reduction on ", temp_word)
            new_temp_word = gog.cascade(temp_word, file_dir, kbmag_ftn_dir)
            #print("after cascading, the temp word is ", temp_word)
            IsCascaded = True
            #if the word changed, free reduction needs to be performed again
            if new_temp_word != temp_word:
                IsFreeReduced = False
        counter +=1
        temp_word = new_temp_word
    #no more reductions can be completed, so return temp_word
    if wd_len(temp_word)==0:
        temp_word = "IdWord"
        #print(temp_word)
    return temp_word

def find_M(gog, file_dir, kbmag_ftn_dir):
    """
    For the graph of groups gog, given that many FSA's have been calculated 
    correctly, calculate the maximum length of a synchronous word difference
    and call this maximum length M. 

    Parameters
    ----------
    gog : a graph of groups object
    
    Variables
    ---------
    M: fellow traveling constant for pi1 of GOG
    K: max(fellow traveling constant for (Gv, Ge) and (Gv, G_rev{e}))
    k2: fellow traveling constant for Ge
    n1: maximum of (number of states of F)+1  and we consider all F automata of the automatic structure of Ge
    
    Returns
    -------
    M : an integer

    """
    Ge_wa = FSA.create_fsa_from_file(file_dir + gog.vertices[0].subgp+".wa")
    Ge_gm = FSA.create_fsa_from_file(file_dir + gog.vertices[0].subgp+".gm")
    Ge_gm_list = gm_split(Ge_gm, file_dir, kbmag_ftn_dir)
    #calculating K
    K=2*max(gog.edges[0].fg.coswa.states.size, gog.edges[0].rg.coswa.states.size)-1
    k3 = 2*(K^2 + 1)
    #calculating upper bound for k2 using Theorem 2.3.2 of ECHLPT
    k2 = 2*Ge_wa.states.size - 1
    #calculating upper bound for n1
    n_gm = max(gm.states.size for gm in Ge_gm_list)
    n1 = max(Ge_wa.states.size, n_gm)+ 1
    print('k2', k2, 'K', K, 'n1', n1)
    #calculating M
    M = max(k3, k2)^2 + max(k3, k2)*n1
    return M
    